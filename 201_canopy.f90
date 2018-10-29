!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_lad()
!--------------------------------------------------------------------!
!    initialize the leaf area density profile                                          
!--------------------------------------------------------------------! 
use types, only:rprec
use param
implicit none
! ---
integer :: jz
! ---
if (use_res_canopy) then
    allocate(a_leafz_uv(nz))
    allocate(a_leafz_w(nz))
    a_leafz_uv = 0._rprec
    a_leafz_w  = 0._rprec

    if (ocean_flag) then
        if (use_field_scale) then
            do jz = 1, nz-1
                if ( canopy_type .eq. 1 ) then              ! suspended canopy
                    if ( zuvp(jz) .ge. hgap .and. zuvp(jz) .le. hgap+hc ) then
                        a_leafz_uv(jz) = lad
                    end if

                    if ( zw(jz) .ge. hgap .and. zw(jz) .le. hgap+hc ) then
                        a_leafz_w(jz) = lad
                    end if
                else if ( canopy_type .eq. 2 ) then         ! submerged / emergent canopy
                    if ( zuvp(jz) .ge. (lz_tot-hc) ) then    
                        a_leafz_uv(jz) = lad
                    end if

                    if ( zw(jz) .ge. (lz_tot-hc) ) then
                        a_leafz_w(jz) = lad
                    end if
                else
                    write(*,*) "Canopy type NOT supported!"
                    stop
                end if
            end do
        end if
    else
        if (use_field_scale) then
            do jz = 1, nz-1
                if ( zuvp(jz) .le. hc ) then
                    a_leafz_uv(jz) = lad
                end if

                if ( zw(jz) .le. hc ) then
                    a_leafz_w(jz) = lad
                end if
            end do
        end if
    end if
end if
! ---
end subroutine init_lad
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine canopy_model 
!-----------------------------------------------------------------------
!   pressure gradient and canopy drag
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
implicit none
! ---
real(kind=rprec) :: V_uv, V_w
integer :: jx, jy, jz
! ---

dragx = 0._rprec
dragy = 0._rprec
dragz = 0._rprec

if (use_res_canopy) then
    if (use_field_scale) then
        do jz = 1, nz - 1
        do jy = 1, nynpy
        do jx = 1, nxt
            V_uv = (  u(jx,jy,jz)**2._rprec &
                 +    v(jx,jy,jz)**2._rprec &
                 +  ((w(jx,jy,jz)+w(jx,jy,jz+1))/2._rprec)**2._rprec)**0.5_rprec
            V_w  = (((u(jx,jy,jz)+u(jx,jy,jz-1))/2._rprec)**2._rprec &
                 +  ((v(jx,jy,jz)+v(jx,jy,jz-1))/2._rprec)**2._rprec &
                 +    w(jx,jy,jz)**2._rprec                         )**0.5_rprec

            dragx(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_uv(jz)*V_uv**(vel_exp-1.)*u(jx, jy, jz) 
            dragy(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_uv(jz)*V_uv**(vel_exp-1.)*v(jx, jy, jz)
            dragz(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_w(jz)*V_w**(vel_exp-1.)*w(jx, jy, jz)
        end do
        end do
        end do

        ! if (rank .eq. 0) then 
        !     open(9,file='../output/drag_check.dat')
        !     write(9,*) 'variables = "x","y","z","dx","dy","dz"'
        !     write(9,*) 'zone i=', nxt, 'j=', nynpy, 'k=', nz, 'f=point'
        !     do jz = 1, nz
        !     do jy = 1, nynpy
        !     do jx = 1, nxt              
        !         write(9, '(6e15.4)') (jx-1) * lx_tot * z_i / nxt, (jy-1) * ly_tot * z_i / nyt, (jz-1) * z_i / (nzt-1), &
        !                             dragx(jx, jy, jz), dragy(jx, jy, jz), dragz(jx, jy, jz) 
        !     end do
        !     end do
        !     end do
        !     close(9)    
        ! end if
    end if 
end if
! ---
end subroutine canopy_model
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------

    



