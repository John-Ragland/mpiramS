subroutine ram(zsrc,rg)

! frq		frequency
! zsrc		source depth
! rg		vector of output ranges
! dr    	range step
! zmax  	max depth
! dz    	depth grid increment
! c0    	"mean" sound speed
! np    	# of pade coefficients
! ns    	# of stability terms
! rs    	stability range
! rb		bathymetry range 
! zb		bathymetry
! rp 		profile ranges(nrg)
! zw    	sound speed grid depth(nzw)
! cw    	sound speed(nzw,nrg)
! zs    	sediment speed grid depth(nzs)
! cs    	sediment speed(nzs,nrg)
! zr		density depth grid(nzr)
! rho		density(nr,nzr)
! za		attenuation depth grid(nza)
! attn		attenuation(nr,nza)

use kinds
use interpolators
use envdata
use mattri
use param
use profiles
use fld

implicit none

real(kind=wp),dimension(:),intent(in) :: zsrc,rg 

integer :: iflag
integer :: ii,n,nz,nr,nb,ir,irl,izl,izll,irr,upd,izs
integer :: ir0(1)
real(kind=wp) :: omega, dr, rend, rnow, rint, rsc
real(kind=wp) :: delzs,zbc,zsc
real(kind=wp) :: rint1(1),zbc1(1), maxrb1

real(kind=wp),dimension(:),allocatable :: rb1, zb1

  omega=2.0_wp*pi*frqq
  k0=omega/c0

! cfact, a1,a2,a3 are used in matrc.  One need not recompute them every time!
! These are double precision because matrc seems to need that for small deltaz.
  cfact=0.5_wp2/deltaz**2
  a1=k0*k0/6.0_wp2
  a2=2.0_wp2*k0*k0/3.0_wp2
  a3=a1

  nz=floor(zmax/deltaz-0.5_wp)
  if (.not.allocated(zg)) call linspace(zg, 0.0_wp, zmax, nz+2)

  nr=size(rg)
  rend=rg(nr)
  rnow=0.0_wp
  if (nr>1) rnow=rg(1)

  dr=deltar    ! dr is deltar of peramx; 
            ! it may need to adjust to get to the range rend precisely
  !dr=abs(dr)
  if ((rend-rnow)<0) dr=-dr

  rsc=abs(rend-rnow)-rs
  if (rs<abs(dr)) rsc=0.0_wp

  ! The initial profiles.
  call profl(rnow+dr/2.0_wp,omega,3)

  ir0=minloc(abs(rp-(rnow+dr/2.0_wp)))
  ir=ir0(1)
  irl=ir

  nb=size(rb)
  allocate(rb1(nb+1),zb1(nb+1))
  ! The initial depth.
  forall(ii=1:nb) rb1(ii)=rb(ii)
  forall(ii=1:nb) zb1(ii)=zb(ii)
  rb1(nb+1)=2.0_wp*rb(nb)+dr
  maxrb1=maxval(rb1)
  zb1(nb+1)=zb(nb)
  rint=rnow+dr/2.0_wp
  rint1(1)=rint      ! for the interpolation, everything must be an array...sigh.
                     ! scalars have rank 0, vectors rank 1
  zbc1=interp1(rb1,zb1,rint1,zb(1))
  zbc=zbc1(1)        ! convert the array zbc1 back to scalar zbc.
  if (rint>maxrb1) zbc=zb1(nb+1)
  iz=floor(1.0_wp+zbc/deltaz); iz=max(2,iz); iz=min(nz,iz)


  ! Self starter
  if (size(zsrc)==1) then
    allocate(uu(nz+2))
    ! zero uu
    uu=0.0_wp*uu
    ! Conditions for the delta function.
    zsc=1.0_wp+zsrc(1)/deltaz
    izs=floor(zsc)
    delzs=zsc-real(izs,wp)
    uu(izs)=(1.0_wp-delzs)*sqrt(2.0_wp*pi/k0)/(deltaz*alpw(izs))
    uu(izs+1)=      delzs* sqrt(2.0_wp*pi/k0)/(deltaz*alpw(izs))

    ! Divide the delta function by (1-X)**2 to get a smooth rhs.
    allocate(pdu(np),pdl(np)); pdu=0.0_wp*pdu; pdl=0.0_wp*pdl
    pdu(1)=cmplx(0.0_wp,wp); pdl(1)=cmplx(-1.0_wp,wp)
    call matrc
    call solvetri
    call solvetri
    deallocate(pdu,pdl)

    ! Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*dr*sqrt(1+X)). (3-D)
    dre=abs(dr)
    ip=0
    call epade
    ! Apply the operator (1-X)**2*(1+X)**(-1/2)*exp(ci*k0*dr*sqrt(1+X)). (2-D)
    call matrc
    call solvetri
    rnow=rnow+dr

    ! reverse propagation backwards to rnow
    dre=-abs(dr)
    ip=1
    call epade
    call matrc
    call solvetri
    rnow=abs(rnow-dr)

    !% lossy propagation backwards to rnow
    ![pdu, pdl]=epade( np, ns, 1, k0, abs(dr));
    ![r1,r2,r3,s1,s2,s3,f3]=matrc(k0,dz,izb,rhob,alpw,alpb,ksqw,ksqb,pdu,pdl);
    !uu=conj(uu);
    !uu=conj(solvetri( izb,uu,r1,r2,r3,s1,s2,s3));
    !rnow=rnow-dr;
  end if

  ! The propagation matrices.
  dre=abs(dr)
  ip=1
  call epade
  call matrc

  !specified starting field
  if (size(zsrc)==size(zg)) uu=zsrc/f3

  if (.not.allocated(psi)) allocate(psi(nz+2,nr))
    psi=0.0_wp*psi
  if (.not.allocated(rout)) allocate(rout(nr))
    rout=0.0_wp*rout

  izll=iz    ! izll is used as a flag for updating profiles for changing bathymetry
! March the acoustic field out in range, inching along at deltar increments.
  do irr=1,size(rg)
    rend=rg(irr)
  
    do 
      if (abs(rnow-rend)<0.1_wp) exit    ! If we're within 10 cm, call it a day.
                                         ! Avoid rounding issues...
      upd=0
      iflag=0

      if (abs(rend-rnow)<abs(dr)) then
        dr=rend-rnow
        dre=abs(dr)
        ip=1
        call epade
        upd=1
      end if

      rint=rnow+dr/2.0_wp
      rint1(1)=rint   ! for the interpolation, everything must be an array...sigh.

      ! Varying bathymetry.
      zbc1=interp1(rb1,zb1,rint1,zb1(1))
      zbc=zbc1(1)     ! convert array zbc1 back to scalar zbc.
      if (rint>maxrb1) zbc=zb1(nb+1)
      izl=iz; iz=floor(1.0_wp+zbc/deltaz); iz=max(2,iz); iz=min(nz,iz)
      if (iz/=izl)  then
      ! bathymetry has changed; call matrc 
          upd=1    
          if (abs(izll-iz)*deltaz > 20.0_wp) then 
          ! The depth has changed by more than 20 m; update the bottom profiles
          ! This is mainly for attenuation and density.
          ! Don't need to call this for EVERY depth change! (I don't think...)
             iflag=iflag+1
             izll=iz
          end if
      end if

      ! Varying profiles - using profile closest to present range.
      ir0=minloc(abs(rp-rint))
      irl=ir; ir=ir0(1)
      if (ir/=irl) then 
      ! sound speed has changed; update profiles and call matrc
           iflag=iflag+2
           upd=1    
      end if

      ! Turn off the stability constraints if the range left is less than rsc.
      if (abs(rend-rnow)<rsc) then
        ns=0   ! ns=1
        rsc=0.0_wp
        upd=1
        dre=abs(dr)
        ip=1
        call epade
      end if

      if (upd==1) then
         if (iflag/=0) call profl(rint,omega,iflag)  ! iflag=1, update bottom properties
         call matrc                                  ! iflag=2, update water sound speed
      end if                                         ! iflag=3, update both

      call solvetri

      rnow=rnow+dr;
    end do

    psi(:,irr)=uu*f3
    rout(irr)=rnow
  end do

    deallocate(uu,r1,r2,r3,s1,s2,s3,f3)
    deallocate(cwg,rhob,alpw,alpb,ksqw,ksqb,pdu,pdl)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine profl(r,omega,iflag)

use kinds
use envdata
use profiles

implicit none

real(kind=wp), parameter :: eta=0.018323389971986   !  eta=1/(40*pi*log10(exp(1)));
                                 !  Why calculate a division, a log, and an exponential every time?
complex(kind=wp), parameter :: ci=cmplx(0.0, 1.0,wp)

integer, intent(in) :: iflag  ! iflag=3 update all; iflag=1 update bathymetry; iflag=2 update sound speed
real(kind=wp), intent(in) :: r, omega

integer :: ir,ii
integer :: ir0(1)
real(kind=wp) :: depth
real(kind=wp), dimension(:) :: rwork(1),zwork(4)
real(kind=wp), dimension(:,:) :: work(1,1)
real(kind=wp), dimension(:), allocatable :: csg,attng

  n=size(zg)
  if (.not.allocated(rhob)) allocate(cwg(n),rhob(n),ksqw(n),alpw(n),alpb(n),ksqb(n))
         ! These variables are stored in the module profiles.f90

if (iflag==2.or.iflag==3) then    ! update water sound speed
    ir0=minloc(abs(rp-r)); ir=ir0(1)
    cwg=gorp2(zw,cw(:,ir),zg)     ! gorp2 uses cubic spline to interpolate to deltaz grid
                                  ! The range dependence is to just use the profile nearest to r.
end if

if (iflag==1.or.iflag==3) then    ! update sediment sound speed, density, attenuation
    allocate(csg(n),attng(n))

!   First find the depth at this range
    rwork(1)=r
    work(:,1)=interp1(rb,zb,rwork,zb(1))
    depth=work(1,1)
    ! The four values of depth that go with cs, rho, and attn
    zwork(1)=0.0_wp; zwork(2)=depth
    zwork(3)=depth+sedlayer; zwork(4)=max(zg(n),zwork(3)+1.0E-6)

! Set up sediment sound speed to increase linearly below the sea floor, with
! a sedlayer-m thick sediment layer. 
    csg=gorp(zwork,cs,zg)
    csg=cwg+csg

! Set up the sediment density and attenuation profiles.
! Attenuation and density follow the bottom.
    rwork(1)=rho(1) 
    rhob=gorp(zwork,rwork,zg)   ! send it only one value, so gorp does the easy thing; rhob a constant.
    attng=gorp(zwork,attn,zg)
end if

! ic0 is 1/c0, calculated in peramx; it doesn't change.
  if (iflag==2.or.iflag==3) then
    forall(ii=1:n) ksqw(ii)=(omega/cwg(ii))**2-(omega*ic0)**2
    forall(ii=1:n) alpw(ii)=sqrt(cwg(ii)*ic0)
  end if
  if (iflag==1.or.iflag==3) then
    forall(ii=1:n) ksqb(ii)=((omega/csg(ii))*(1.0_wp+ci*eta*attng(ii)))**2-(omega*ic0)**2
    forall(ii=1:n) alpb(ii)=sqrt(rhob(ii)*csg(ii)*ic0)
    deallocate(csg,attng)
  end if

  !csg=cwg.*sqrt(1+(eta*attng).^2)./(1+ci*eta*attng);
  !ksqb=((omega./csg)).^2-(omega/c0)^2;
  !alpb=sqrt(rhob.*abs(csg)/c0);
end subroutine profl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gorp(x,y,u)
! Sets up a simple linear interpolation

  use kinds

  implicit none

  real(kind=wp),dimension(:),allocatable :: gorp
  real(kind=wp),dimension(:) :: x,y,u

  integer :: n,m,ii

  n=size(u)
  m=size(y)   ! x and y the same size

  allocate(gorp(n))
  
  if (size(y)==1) then
    !forall(ii=1:n) gorp(ii)=y(1)
    gorp=y(1)+0.0_wp*gorp
    return
  end if

  gorp=interp1(x,y,u,y(1))
  
  if (x(m)>=u(n)) return

  do ii=1,n
    ! assume monotonically increasing x
    if (u(ii)>x(m)) gorp(ii)=y(m)
  end do

end function gorp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gorp2(x,y,u)
! Cubic spline interpolation using cubspl/ppvalu

  use kinds

  implicit none

  real(kind=wp),dimension(:) :: x,y,u

  real(kind=wp),dimension(:),allocatable :: gorp2
  real(kind=wp2),dimension(:,:),allocatable :: cry
  real(kind=wp2) :: ztemp,ctemp

  integer :: n,m,ii

interface
  subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
    use kinds
    IMPLICIT NONE
    integer :: ibcbeg,ibcend,n
    real ( kind = wp2 ) c(4,n)
    real ( kind = wp2 ) tau(n)
  end subroutine cubspl

  function ppvalu (break, coef, l, k, x, jderiv )
    use kinds
    IMPLICIT NONE
    integer :: jderiv,k,l
    real ( kind = wp2 ) ppvalu 
    real ( kind = wp2 ) break(l+1)
    real ( kind = wp2 ) coef(k,l)
    real ( kind = wp2 ) x
  end function ppvalu
end interface

  n=size(u)
  m=size(y)   ! x and y the same size

  allocate(gorp2(n))
  
  if (size(y)==1) then
    forall(ii=1:n) gorp2(ii)=y(1)
    return
  end if

  !gorp2=interp1(x,y,u,y(1))
  allocate(cry(4,m))
  cry(1,:)=real(y,wp2)
  call cubspl(real(x,wp2),cry,m,0,0)
! Returns coefficients for obtaining sound speed, dc/dz and d2c/dz2.
! These things are used by the ppvalu function.
! Evaluate sound speed at all depths.
  do ii=1,n
     ztemp=u(ii)
     ctemp=ppvalu(real(x,wp2),cry,m-1,4,ztemp,0)
! According to gprof, 3% of the calculation time is in this call to ppvalu.
     gorp2(ii)=real(ctemp,wp)
  enddo

  if (x(m)>=u(n)) return

  do ii=1,n
    ! assume monotonically increasing x
    if (u(ii)>x(m)) gorp2(ii)=y(m)
  end do

end function gorp2

end subroutine ram
