!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine project()
!-----------------------------------------------------------------------
!--provides u, v, w at 1:nz
!-----------------------------------------------------------------------
use param
use sim_param
! use scalars_module, only:inflow_pollen
implicit none
! ---
integer :: jx, jy, jz, jz_min
real(rprec) :: RHS
! ---
do jz = 1, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    RHS = -tadv1 * dpdx(jx, jy, jz)
    u(jx, jy, jz) = u(jx, jy, jz) + dt * RHS
    RHS = -tadv1 * dpdy(jx, jy, jz)
    v(jx, jy, jz) = v(jx, jy, jz) + dt * RHS
end do
end do
end do

if ( coordz == 0 ) then
    jz_min = 2
else
    jz_min = 1
end if

do jz = jz_min, nz - 1
do jy = 1, nynpy
do jx = 1, nxt
    RHS = -tadv1 * dpdz(jx, jy, jz)
    w(jx, jy, jz) = w(jx, jy, jz) + dt * RHS
end do
end do
end do

! if (inflow) call inflow_cond
! ! Zero inflow condition for pollen Chamecki 08/21/2006
! if (inflow_PCon) call inflow_pollen()

!--enfore bc at top
if ( coordz == npz - 1 ) then
    w(:, :, nz) = 0._rprec
end if

if ( coordz == 0 ) then
    w(:, :, 1) = 0._rprec
end if
! ---
end subroutine project
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! subroutine inflow_cond()
! !-----------------------------------------------------------------------
! !
! !-----------------------------------------------------------------------    
! use types, only:rprec
! use param, only:face_avg, jx_s, nxt, nyt, nz, pi, read_inflow_file, &
!                 passive_scalar, x_relax
! use sim_param, only:u, v, w, theta
! use io, only:inflow_read
! implicit none
! ! ---
! integer :: i, jx, jx_0, jx_s1
! real(rprec) :: factor

! !--read from file
! if (read_inflow_file) then !--read vel inflow @ jx = jx_s from file
!     call inflow_read() !--this sets u, v, w at (jx_s,:,:)
! else
!     u(jx_s, :, :) = face_avg
!     v(jx_s, :, :) = 0._rprec
!     w(jx_s, :, :) = 0._rprec
! end if

! !--just an experiment
! jx_s1 = modulo(jx_s + 1 - 1, nxt) + 1
! u(jx_s1, :, :) = u(jx_s, :, :)
! v(jx_s1, :, :) = v(jx_s, :, :)
! w(jx_s1, :, :) = w(jx_s, :, :)

! jx_0 = jx_s - x_relax !--velocity at this point is not forced here

! do i = 1, x_relax - 1

!     jx = jx_s - i

!     factor = 0.5_rprec*(1._rprec - cos(pi*real(i, rprec)/x_relax))
!     !factor = real (i, rprec) / x_relax

!     u(jx, 1:nyt, 1:nz) = u(jx_s, 1:nyt, 1:nz) + &
!         factor*(u(jx_0, 1:nyt, 1:nz) - u(jx_s, 1:nyt, 1:nz))
!     v(jx, 1:nyt, 1:nz) = v(jx_s, 1:nyt, 1:nz) + &
!         factor*(v(jx_0, 1:nyt, 1:nz) - v(jx_s, 1:nyt, 1:nz))
!     w(jx, 1:nyt, 1:nz) = w(jx_s, 1:nyt, 1:nz) + &
!         factor*(w(jx_0, 1:nyt, 1:nz) - w(jx_s, 1:nyt, 1:nz))

!     if (passive_scalar) then
!         theta(jx, 1:nyt, 1:nz) = 0._rprec + factor*(theta(jx_0, 1:nyt, 1:nz) - 0._rprec)
!     end if

! end do

! end subroutine inflow_cond
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! subroutine inflow_pollen()
! ! This is a simplified version of inflow_cond() in forcing.f90
! ! Chamecki 08/21/2006
! use types, only : rprec
! use param, only : jx_p, nyt, nz, pi, x_relaxp
! use sim_param, only : pcon
! implicit none

! integer :: jx
! integer :: i
! integer :: jx_0
! integer :: ipcon

! real (rprec) :: factor

! ! Initial point for damping
! jx_0 = jx_p - x_relaxp

! ! Changed from i=1 to i=0 to include last point in the domain
! do i = 0, x_relaxp-1
!     jx = jx_p - i
!     factor = 0.5_rprec * (1._rprec - cos (pi * real (i, rprec) / x_relaxp))
!     do ipcon=1,npcon
!         PCon(jx, 1:nyt, 1:nz, ipcon) =  factor * PCon(jx_0, 1:nyt, 1:nz, ipcon)
!     end do
! end do
! ! ---
! end subroutine inflow_pollen
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
