module topbc
!--------------------------------------------------------------------!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!     
use types, only: rprec
use param, only: nz, nzt, damping_method
implicit none
save
! ---
real(kind=rprec), dimension(:), allocatable :: sponge

contains

    subroutine setsponge()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!     
    use param
    implicit none
    ! ---
    real(kind=rprec) :: factor
    real(kind=rprec) :: z_d, cfrdmp, z_local
    integer :: k
    
    integer :: k_global
    real(rprec) :: sponge_top
    ! ---
    if ( damping_method == 2 ) then
        z_d = 0.75_rprec * lz_tot
        cfrdmp = 3.9_rprec
        sponge = 0._rprec
        do k = 1, nz - 1
            if (zuvp(k) .ge. z_d .and. zuvp(k) .le. lz_tot) then
                sponge(k) = cfrdmp*0.5_rprec*(1._rprec - cos(pi*(zuvp(k) - z_d)/(lz_tot - z_d)))
            else
                sponge(k) = 0._rprec
            end if
        end do

        if ( coordz == npz - 1 ) then
            sponge(nz) = sponge(nz - 1)
        end if
    else if (damping_method == 1) then
    ! sets relaxation term to vertical momentum equation in top quarter
    ! of domain, relaxation timescale on top is 50s with a factor of 5 if we
    ! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear
    ! flows
        sponge = 0._rprec
        factor = 9._rprec / (nzt - 3*nzt/4 + 1)
        sponge_top = z_i / (50._rprec * u_scale)
        do k = 1, nz-1
            k_global = k + coordz * (nz - 1)
            if (k_global > 3*nzt/4 + 1) then
                sponge(k) = sponge_top * 5._rprec**((k_global-nzt) * factor)
            end if
        end do
    end if
    ! ---
    end subroutine setsponge
! ---    
end module topbc
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine add_damping_layer()
!--------------------------------------------------------------------!
!   Add damping terms to the momentum RHS                                  
!--------------------------------------------------------------------!     
use param, only: nxt, nynpy, nz, ubc_mom, ubc, damping_method
use sim_param, only: u, v, w, RHSx, RHSy, RHSz
use topbc, only: sponge 
implicit none
! ---
integer :: jz
! --- 
if ( ubc_mom == 'stress free' ) then
    if (ubc == 1 .and. damping_method == 2) then
        do jz = 1, nz - 1
            RHSx(1:nxt, 1:nynpy, jz) = RHSx(1:nxt, 1:nynpy, jz) - 0.5*(sponge(jz) + sponge(jz + 1))* &
                                        (u(1:nxt, 1:nynpy, jz) - sum(u(1:nxt, 1:nynpy, jz))/(nxt*nynpy))
            RHSy(1:nxt, 1:nynpy, jz) = RHSy(1:nxt, 1:nynpy, jz) - 0.5*(sponge(jz) + sponge(jz + 1))* &
                                        (v(1:nxt, 1:nynpy, jz) - sum(v(1:nxt, 1:nynpy, jz))/(nxt*nynpy))
            RHSz(1:nxt, 1:nynpy, jz) = RHSz(1:nxt, 1:nynpy, jz) - 0.5*sponge(jz)* &
                                        (w(1:nxt, 1:nynpy, jz) - sum(w(1:nxt, 1:nynpy, jz))/(nxt*nynpy))
        end do
    else if (ubc == 1 .and. damping_method == 1) then
        do jz = 1, nz - 1
            RHSz(1:nxt, 1:nynpy, jz) = RHSz(1:nxt, 1:nynpy, jz) - sponge(jz) * w(1:nxt, 1:nynpy, jz)
        end do
    end if
end if
! ---
end subroutine add_damping_layer
