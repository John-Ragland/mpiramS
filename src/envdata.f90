module envdata
! sound speed and bathymetry data and parameters
! Modules such as these take the place of "common" blocks, but
! are meant to be better - less error prone.

use kinds

implicit none

real(kind=wp) :: sedlayer
real(kind=wp),dimension(:),allocatable :: zg               ! depth grid with deltaz spacing
real(kind=wp),dimension(:),allocatable :: rb,zb            ! bathymetry
real(kind=wp),dimension(:),allocatable :: rp,zw            ! range and depths of sound speeds
real(kind=wp),dimension(:,:),allocatable ::  cw            ! sound speeds, etc.
real(kind=wp),dimension(:), allocatable :: cs,rho,attn    
                           ! bottom properties are simple and range independent - four values.

end module envdata

