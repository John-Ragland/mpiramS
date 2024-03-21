program peramx
!
! RAM PE example - now in Fortran 95!
!
! This code was copied wholesale from Matt Dzieciuch's (@SIO) matlab version, and 
! adapted to Fortran 95.  Matt, in turn, developed his code based on Mike Collins' 
! original FORTRAN version.
!
! B. Dushaw  May 2008
! Single precision version July 2013

use kinds
use interpolators
use envdata
use param
use fld

implicit none

real(kind=wp) :: cmin
real(kind=wp),dimension(:),allocatable :: rp0

real(kind=wp),dimension(:),allocatable :: zg1
complex(kind=wp),dimension(:),allocatable :: psi1
complex(kind=wp),dimension(:,:),allocatable :: psif

! input parameters - c.f., file "in.pe"
integer :: dzm, iflat, ihorz, ibot
real(kind=wp) :: fc,Q,T,dum
real(kind=wp),dimension(:),allocatable :: zsrc,rmax
character(len=20) :: name1,name2     ! sound speed and bathymetry filenames
real(kind=wp),dimension(:),allocatable :: eps
real(kind=wp),dimension(:,:),allocatable :: cq

integer :: nss

integer :: nb,nzp,nrp,nrp0,n,nf1,nf
real(kind=wp) :: bw, fs, Nsam, df, tmp
real(kind=wp),dimension(:),allocatable :: frq

integer ::  nzo,icount
real(kind=wp) :: omega

real(kind=wp) :: rate
integer :: t1,t2,cr,cm

integer :: ii,jj,iff,length

integer, parameter :: nunit=2
complex(kind=wp), parameter :: j=cmplx(0.0_wp,1.0_wp)
complex(kind=wp) :: scl

interface
  subroutine ram(zsrc,rg)

    use kinds
    use interpolators
    use envdata
    use fld
    use param

    implicit none

    real(kind=wp),dimension(:),intent(in) :: zsrc,rg

  end subroutine ram

end interface

allocate(zsrc(1),rmax(1))

open(nunit,file='in.pe',status='old')
read (2,'(f4.0)')  fc                ! skip the first line - read a dummy
read (2,'(f4.0,1x,f2.0)') fc,Q       ! center frequency and Q
read (2,'(f4.1)') T              ! time window width
read (2,'(f6.1)') zsrc(1)        ! source depth
read (2,'(f12.3)') rmax(1)       ! receiver range
read (2,'(f5.2)') deltaz         ! depth accuracy parameter
read (2,'(f6.2)') deltar         ! range accuracy parameter
read (2,'(i1,1x,i1)') np,nss     ! np -# pade coefficients
                                 ! ns -# stability terms
read (2,'(f7.1)') rs             ! stability range
read (2,'(i2)') dzm              ! dzm - depth decimation
read (2,'(a20)') name1           ! sound speed filename; "munk" will just use canonical
name1=trim(name1) ! remove trailing blanks
read (2,'(i1)') iflat            ! 0=no flat earth transform, 1=yes
read (2,'(i1)') ihorz            ! 0=no horizontal linear interpolation, 1=yes
read (2,'(i1)') ibot             ! 0=no bottom, 1=bottom and read a file
read (2,'(a20)') name2           ! bathymetry filename; ignored if ibot=0
name2=trim(name2) ! remove trailing blanks

close(nunit)

print *
print '(a)','INPUT PARAMETERS:'
print '(a,f10.2)','Center frequency (Hz): ', fc
print '(a,f2.0)','Q: ', Q
print '(a,f5.2)','Bandwidth (f0/Q - Hz): ', fc/Q
print '(a,f4.1)','Time window width (s): ', T
print '(a,f6.1)','Source depth (m): ', zsrc(1)
print '(a,f12.3)','Range (m): ', rmax(1)
print '(a,f5.2)','Depth accuracy (deltaz, m): ', deltaz
print '(a,f6.2)','Range accuracy (deltar, m): ', deltar
print '(a,i2)','No. of pade coefficients: ', np
print '(a,i2)','No. of stability terms: ', nss
print '(a,f7.1)','Stability range (m): ', rs
print '(a,i2)','Output depth decimation: ', dzm
print '(a,a)','Sound speed filename: ', name1
print '(a,i1)','Flat-earth transform flag: ', iflat
print '(a,i1)','Horizontal interpolation flag: ', ihorz
print '(a,i1)','Ocean bottom flag: ', ibot
print '(a,a)','Ocean bottom filename: ', name2
print *

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SOUND SPEED SET UP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (name1=='munk') then
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! *** set up range independent Munk SSP ***
   nrp=1
   allocate(rp(nrp))   ! must be a vector, since ram expects a vector
   rp(nrp)=0.0_wp

   nzp=5001
   allocate(zw(nzp),cw(nzp,nrp))
   forall (ii=1:nzp) zw(ii)=real(ii-1,wp)
   cw(:,1)=cssprofile(zw,0.0_wp)     ! 2nd parameter is latitude
else
   ! open and read the sound speed from file name1
   print *,'Reading sound speed file.'
   open(nunit,file=name1,status='old')
   ! first run through and count the -1's to get the number of profiles and depths
     nzp=0
     nrp0=0
     do
        read(nunit,*,end=10) dum
        if (dum<0) then
            nrp0=nrp0+1 
        end if
        if ((dum>=0).and.(nrp0==1)) nzp=nzp+1
     end do
 10  print *,'Found ',nrp0,' profiles, and ',nzp,' depths.' 
   close(nunit)
     allocate(rp0(nrp0),zw(nzp),cq(nzp,nrp0))
   open(nunit,file=name1,status='old')
     do ii=1,nrp0
        read(nunit,*) dum,rp0(ii)
        rp0(ii)=rp0(ii)*1000.0_wp
        do jj=1,nzp
            read(nunit,*) zw(jj),cq(jj,ii)
        enddo
     end do 
   close(nunit)

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Sound speeds are read in. 
   ! Now interpolate in range to 10-km intervals (linear), say, if requested.
   
   if (ihorz==1) then
   ! Horizontal
     nrp=nint(rmax(1)/10000.0_wp)
     call linspace(rp, rp0(1),rmax(1),nrp)
     allocate(cw(nzp,nrp))
     do jj=1,nzp
        cw(jj,:)=interp1(rp0,cq(jj,:),rp,cq(jj,1))
     enddo
   else
     nrp=nrp0
     allocate(rp(nrp),cw(nzp,nrp))
     rp=rp0
     cw=cq
   end if
   deallocate(cq,rp0)

     ! flat earth transformation
   if (iflat==1) then
      print *
      print *,'Applying flat-earth transformation.'
      print *
      allocate(eps(nzp))
      eps=zw*invRe
      zw=zw*(1.0_wp+(1.0_wp/2.0_wp)*eps+(1.0_wp/3.0_wp)*eps*eps)
      forall(ii=1:nrp) cw(:,ii)=cw(:,ii)*(1+eps+eps*eps)
      deallocate(eps)
      ! also transform the source depth
      zsrc(1)=zsrc(1)*(1.0_wp+(1.0_wp/2.0_wp)*zsrc(1)*invRe +   &
            (1.0_wp/3.0_wp)*(zsrc(1)*invRe)**2)
   end if

endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! mean sound speed
n=size(cw)
c0=sum(cw)/n
ic0=1.0_wp/c0
cmin=minval(cw)     ! minimum sound speed for calculating tdelay

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOTTOM SET UP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (ibot==0) then
   ! use default flat bottom
   nb=1
   allocate(rb(nb),zb(nb))
   rb(nb)=0.0_wp
   zb(nb)=maxval(zw) - 400.0_wp    
    ! Set the ocean floor to be 400 m shallower than deepest sound speed depth.
    ! c.f., attenuation below
else
   ! open and read in the bathymetry from file name2
   print *,'Reading bathymetry file.'
   open(nunit,file=name2,status='old')
! first find the number of values
   nb=0
   do
     read(nunit,*,end=15) dum
     nb=nb+1 
   end do
15 print *,'Found ',nb,' bathymetry values.' 
   rewind(nunit)
   allocate(rb(nb),zb(nb))
   do ii=1,nb
      read(nunit,*) rb(ii),zb(ii)
   end do 
   close(nunit)
endif

   ! flat earth transformation
if (iflat==1) then
   allocate(eps(nb))
   eps=zb*invRe
   zb=zb*(1.0_wp+(1.0_wp/2.0_wp)*eps+(1.0_wp/3.0_wp)*eps*eps)
   deallocate(eps)
   !do ii=1,nb
   !   if (zb(ii)>maxval(zw)) then
   !      print *,'Bottom deeper than water sound speed',zb(ii),maxval(zw)
   !      zb(ii)=maxval(zw)
   !   end if
   !end do  
end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOTTOM PROPERTIES MODEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Sediment speed, density and attenuation are gross approximations.
!  (zw,cw) is depth, water sound speed
!   and (zs,cs) is depth, sediment sound speed

! For the bottom we introduce a simple range-indendent model, but one which
! follows the bathymetry.  The profl subroutine in ram.f90 implements this;
! here the basic model is given.
! Bottom speed, density, and attenuation follow rb,zb

! Sediment layer thickness - meters
sedlayer=300.0_wp

! Sediment sound speed - this will be speed relative to the water sound speed.
! Four values are given:  at the surface, at the bottom, sedlayer-m below the bottom,
! and at the center of the earth...
allocate(cs(4))
cs(1)=0.0_wp   
cs(2)=0.0_wp
cs(3)=200.0_wp 
cs(4)=200.0_wp  ! sound speed increases by 200 m/s sedlayer-m
                                 ! below the sea floor, that is to about 1700 m/s.
! One doesn't want an abrupt change in sound speed, or you get too much reflection.
! The sharper the sound speed gradient, the greater the reflection at the sea floor.

! Sediment density - use only one approximate value:  1.2 seems to be nominal.
allocate(rho(4))
rho(1)=1.2_wp 
rho(2)=1.2_wp
rho(3)=1.2_wp 
rho(4)=1.2_wp

! Sediment attenuation - bottom attenuation is 0.5 dB/wavelength from the ocean 
! surface to the sea floor and increases to 5.0 dB/wavelength sedlayer-m below 
! the ocean floor.   0.5 dB/wavelength is a very large value for low frequency.
allocate(attn(4))
attn(1)=0.5_wp  
attn(2)=0.5_wp
attn(3)=5.0_wp  
attn(4)=5.0_wp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                Example using fc=75 Hz, Q=2, and time window T=2.0 
bw=fc/Q        ! 75/2=37.5 Hz bandwidth
fs=4.0_wp*fc   ! 4*75 Hz = 300 Hz sampling frequency
!dt=1.0_wp/fs   ! 1/300 s sampling interval
Nsam=fs*T      ! 300 Hz*2 s = 600 samples
df=fs/Nsam     ! 300/600 = 0.5 Hz frequency interval

! frequency set up
tmp=(bw-df)/df
nf1=int(tmp) + 1  ! (37.5-0.5)/0.5 + 1 = 75 frequencies, one half
nf=2*nf1+1          ! including 0, 151 frequencies altogether.
allocate(frq(nf))
do ii=1,nf1
  frq(ii)=-(nf1-(ii-1))*df + fc
enddo
frq(nf1+1)= fc
do ii=1,nf1
  frq(ii+nf1+1)=ii*df + fc
enddo

! frq=[df:df:bw];
! frq=[-fliplr(frq) 0 frq]+fc;
! nf=size(frq)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zmax=maxval(zw)

! icount is number of depths in zg, psi
icount=floor(zmax/deltaz-0.5_wp)+2

nzo=0
do ii=1,icount,dzm
   nzo=nzo+1 
end do
! nzo is number of depths in psi1, psif, zg1

allocate(psif(nzo,nf))

call system_clock(count_rate=cr)
call system_clock(count_max=cm)
rate=real(cr)

print *,nf,' total frequencies'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Meat and Potatoes
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!$OMP PARALLEL PRIVATE (psi1,omega,scl,t1,t2,cr,rate) 
allocate(psi1(nzo))
!$OMP DO SCHEDULE(STATIC,1)
  do iff=1,nf
    call system_clock(t1)

    frqq=frq(iff)
    ns=nss

    call ram(zsrc,rmax)

! The miracle of fortran95!
    psi1=psi(1:icount:dzm,1)

    omega=2.0_wp*pi*frqq 
    ! 3-D
    scl=exp(j*(omega/c0*rout(1) + pi/4.0_wp))/4.0_wp/pi
    ! 2-D
    ! k0=omega/c0
    !scl=j*exp(j*omega/c0*rout)/sqrt(8.0_wp*pi*k0)
    psif(:,iff)=scl*psi1

    call system_clock(t2,cr)
    rate=real(cr)
    print '(a,i4,1x,f8.3,a,f9.3,a)','iff,frqq = ', iff,frqq, &
            '  Time: ',nint(100.0_wp*dble(t2-t1)/rate)/100.0_wp,' sec'
  end do
!$OMP END DO
!$OMP END PARALLEL

  allocate(zg1(nzo))
  zg1=zg(1:icount:dzm)
!  Remove the flat-earth transform (or most of it, anyways)
  if (iflat==1) then
    allocate(eps(nzo))
    eps=zg1*invRe
    zg1=zg1/(1.0_wp+(1.0_wp/2.0_wp)*eps+(1.0_wp/3.0_wp)*eps*eps)
    deallocate(eps)
  end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! At this point we're done.
! Now need to save the information for plotting in matlab.
! Need to save:  psif(:,:), rout, c0, zg1(:), fs, wind(:), frq(:), Nsam, nf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inquire(iolength=length) real(psif(1,:))
length=2*length+length/nf  ! We need nf real and nf imaginary values and the depth.
! length is the record length, an important piece of information
! Write it out for handy reference:
open(nunit, form='formatted',file='recl.dat')
write(nunit,*) length
close(nunit)

! Now write out the data to a direct access file:
open(nunit, access='direct',recl=length,file='psif.dat')

! Float the integers to real, otherwise there will be trouble.
write(nunit,rec=1) Nsam,real(nf,wp),real(nzo,wp),rout,c0,cmin,fs,Q
write(nunit,rec=2) frq     ! vector of size nf

! Remove the flat-earth transform (or most of it, anyways)
allocate(eps(nzo))
eps=zg1*invRe
zg1=zg1/(1.0_wp+(1.0_wp/2.0_wp)*eps+(1.0_wp/3.0_wp)*eps*eps)
deallocate(eps)

do ii=1,nzo
 write(nunit,rec=ii+2) zg1(ii),((real(psif(ii,jj))),(aimag(psif(ii,jj))),jj=1,nf)
end do

close(nunit)

print *, 'ALL DONE! '
print *,'recl.dat has the record length of the direct access file.'
print *,'The direct access file with the parameters and results is psif.dat.'

stop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cssprofile(z,lat)
! canonical sound speed profile
!
! z (meters)
! c (meters/second)
!
! if abs(lat)<67.5 temperate ocean (default)
! else polar ocean

use kinds

implicit none

real(kind=wp),dimension(:),allocatable :: cssprofile

real(kind=wp) :: lat
real(kind=wp),dimension(:) :: z 
real(kind=wp),dimension(:),allocatable :: eta

integer :: i,n
real(kind=wp) :: gammaa, za, h, ca

gammaa=0.0113_wp/1000.0_wp
za=1000.0_wp
h=1000.0_wp
ca=1500.0_wp

n=size(z)
allocate(cssprofile(n),eta(n))

if (abs(lat)<67.5_wp) then
  !temperate ocean
  do i=1,n
    eta(i)=2.0_wp*(za-z(i))/h
    cssprofile(i)=ca*(1.0_wp+h*gammaa*(exp(eta(i)) - eta(i) - 1.0_wp)/2.0_wp)
  end do
else
  !polar ocean
  do i=1,n
    cssprofile(i)=ca*(1+gammaa*z(i))
  end do
end if

end function cssprofile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end program peramx

