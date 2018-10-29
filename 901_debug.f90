!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
subroutine screen_diagnose
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
use types, only:rprec
use param
implicit none
! ---
integer, parameter :: wbase = 1 ! controls the frequency of screen diagnostics
! ---
if (modulo(jt, wbase) == 0) then
    call check_rmsdiv()             ! calculate the divergence of the flow field
    !call output_loop()
    !  call output_slice_loop()
end if
! ---    
end subroutine screen_diagnose
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!  
subroutine check_rmsdiv()
!--------------------------------------------------------------------!
! actually, this is NOT the rms divergence of velocity, its like an
! l_1 norm or something.
!--------------------------------------------------------------------!  
use types, only: rprec
use param
use sim_param, only: du=>dudx, dv=>dvdy, dw=>dwdz
implicit none
! ---  
integer :: jx, jy, jz, jz_max
real(kind=rprec) :: rms, rms_global
! ---
if ( coordz == npz-1 ) then
    jz_max = nz-1
else
    !jz_max = nz
    jz_max = nz-1
end if

rms = 0._rprec
do jz = 1, jz_max
do jy = 1, nynpy
do jx = 1, nxt
    rms = rms + abs(du(jx,jy,jz)+dv(jx,jy,jz)+dw(jx,jy,jz))
end do
end do
end do
rms = rms / (nxt*nyt*(jz_max))  ! should be nyt here

call mpi_reduce(rms, rms_global, 1, MPI_RPREC, MPI_SUM, 0, comm, ierr)
if (rank == 0) then
    ! rms = rms_global / np
    open(13,file=path//'output/check_div.out',status="unknown",position="append")
    write(13, *) jt, rms
    close(13)
end if
! ---
end subroutine check_rmsdiv
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
subroutine check_cfl(CFL, visc_stab)
!--------------------------------------------------------------------!
! This subroutine computes CFl and viscous stability and is called 
! every wbase timesetps from main.f90
!--------------------------------------------------------------------!
use types, only: rprec
use param
use sim_param, only: u, v, w
use sgsmodule, only: Nu_t
implicit none
! ---
real(kind=rprec), parameter :: CFL_max = 0.6_rprec     ! threshold value of CFL
real(kind=rprec), parameter :: CFL_min = 0.001_rprec
! ---
real(kind=rprec), intent(out) :: CFL, visc_stab
real(kind=rprec) :: u_max, v_max, w_max
real(kind=rprec) :: temp, nu_max
real(kind=rprec) :: cflx, cfly, cflz
real(kind=rprec) :: cflt, cflp, nu_max0

integer :: i, j, k
! ---
cflp = 0._rprec
u_max  = maxval(u(1:nxt, 1:nynpy, 1:nz - 1))
v_max  = maxval(v(1:nxt, 1:nynpy, 1:nz - 1))
w_max  = maxval(w(1:nxt, 1:nynpy, 1:nz - 1))
nu_max = maxval(Nu_t(1:nxt, 1:nynpy, 1:nz - 1))

cflx = u_max * dt / dx
cfly = v_max * dt / dy
cflz = w_max * dt / dz
nu_max0 = nu_max * dt / min(dx, dy, dz)**2

cflt = cflx + cfly + cflz

! do k = 1, nz-1
! do j = 1, nynpy
! do i = 1, nxt
!     cflx = abs(u(i,j,k) * dt / dx)
!     cfly = abs(v(i,j,k) * dt / dy)
!     cflz = abs(w(i,j,k) * dt / dz)
!     cflt = cflx + cfly + cflz
!     cflp = max(cflt, cflp)

!     ! if ( cflp .ge. CFL_max ) then
!     !     write(*,'(i4.4, 6e11.4)') rank, cflx, cfly, cflz, u_max, v_max, w_max
!     ! end if

! end do
! end do
! end do

call mpi_allreduce(cflt, CFL, 1, mpi_rprec, mpi_max, comm, ierr)
call mpi_allreduce(nu_max0, visc_stab, 1, mpi_rprec, mpi_max, comm, ierr)
! ---
if ( CFL .gt. CFL_max ) then
    if ( rank == 0 ) then
        write(*,*) "CFL = ", CFL, "IS TOO LARGE, PLEASE REDUCE THE TIME STEP"
    end if
    call mpi_finalize(ierr)
    stop
!else if ( CFL .lt. CFL_min ) then
!    if ( rank == 0 ) then
!        write(*,*) "CFL = ", CFL, "IS TOO SMALL, PLEASE INCREASE THE TIME STEP"
!    end if
!    call mpi_finalize(ierr)
!    stop
end if

!if ( rank == 0 ) then
!    open(13,file=path//'output/check_cfl.out',status="unknown",position="append")
!    write (13, *) cfl, visc_stab
!    close(13)
!end if
! ---
end subroutine check_cfl
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!  
!subroutine output_loop
!!--------------------------------------------------------------------!
!!BC revised by Bicheng Chen to add output sample
!!-BC add variable flagVSpl, ttVSpl, flagCSpl, ttCSpl
!!-BC add new subroutine checkVSpl and checkCSpl                                               
!!--------------------------------------------------------------------! 
!use param,only:path, output, c_count, theta_flag, theta_init_time, jt, jt_total,  &
!               jan_diurnal_run, flagVSpl, ttVSpl, flagCSpl,             &
!               ttCSpl, flag_srfV, tt_srfV, PCon_FLAG, PCon_init,        &
!               use_avgslice, base_time
!use scalars_module2,only:scalar_slice, pollen_slice, budget_TKE_scalar
!use io, only:calc_mean, avgslice, checkpoint, checkVSpl, checkCSpl,     &
!             check_srfV, post_spec, io_spec, spec_write_start,          &
!             spec_write_end, spec_write_freqz
!implicit none
!! ---
!jt_total = jt_total + 1
!!call calc_mean()
!
!!cyan if (output) then
!!    if (mod(jt_total, base_time)==0) then
!!    !-- move all stuff into checkpoint (Bicheng Chen 06/26/2015)
!!        call checkpoint()
!!    end if
!!end if
!! --- BC added by Bicheng Chen for output the sample data
!if (flagVSpl) then
!    if (mod(jt_total, ttVSpl)==0) then
!        call checkVSpl()
!    end if
!end if
!
!if (flagCSpl) then
!    if (mod(jt_total, ttCSpl)==0) then
!        call checkCSpl()
!    end if
!end if
!!BC END
!! ---
!if (flag_srfV) then
!    if (mod(jt_total, tt_srfV)==0) then
!        call check_srfV()
!    end if
!end if
!
!! ---
!!cyanif ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then
!!!  if ((io_spec) .and. mod(jt_total,spec_write_freqz)==0) then
!!    if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
!!end if
!
!!cyan if (time_spec.gt.0) call timeseries_spec
!! ---
!end subroutine output_loop
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine check_field(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, p, theta, PCon
implicit none
! ---
integer, intent(in) :: num
real, dimension(ldx, nyt, nzt) :: u_tot
real, dimension(ldx, nyt, nzt) :: v_tot
real, dimension(ldx, nyt, nzt) :: w_tot
real, dimension(ldx, nyt, nzt) :: p_tot
real, dimension(ldx, nyt, nzt) :: theta_tot
! real, dimension(ldx, nyt, 1:nzt, npcon) :: PCon_tot

character(80) :: fname
integer :: jx, jy, jz, ipcon
real(rprec) :: xx, yy, zz

! ---
call collect_data_3d(u, u_tot)
call collect_data_3d(v, v_tot)
call collect_data_3d(w, w_tot)
call collect_data_3d(p, p_tot)
call collect_data_3d(theta, theta_tot)    
! call collect_data4_3d(PCon(:,:,:,1), PCon_tot(:,:,:,1))

if ( rank == 0 ) then
    write(fname, '(a,i8.8,a)') '../output/z/vel_1D_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'  
    do jz = 1, nzt-1
       !zz = (jz-1) * z_i / (nzt-1)
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')                      & 
                           sum(u_tot(:, :, jz))/float(nxt*nyt)*u_scale,     & 
                      ',', sum(v_tot(:, :, jz))/float(nxt*nyt)*u_scale,     &
                      ',', sum(w_tot(:, :, jz))/float(nxt*nyt)*u_scale,     &
                      ',', sum(theta_tot(:, :, jz))/float(nxt*nyt)*T_scale
    end do
    close(9)
    
    write(fname, '(a,i8.8,a)') '../output/xz/vel_xz_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    jy = nyt/2
    do jz = 1, nzt-1
    do jx = 1, nxt
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')      &
            u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
            w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale
    end do
    end do
    close(9)
    
    write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_b_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    jz = 24
    do jy = 1, nyt
    do jx = 1, nxt
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
            u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
            w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale 
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_m1_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    jz = 48
    do jy = 1, nyt
    do jx = 1, nxt
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
            u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
            w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale 
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_m2_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    jz = 72
    do jy = 1, nyt
    do jx = 1, nxt
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
            u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
            w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/xy/vel_xy_u_check_', num, '.csv'
    open(9,file=fname)
    write(9, '(7a)') 'u', ',', 'v', ',', 'w', ',', 'theta'
    jz = 96
    do jy = 1, nyt
    do jx = 1, nxt
        write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
            u_tot(jx, jy, jz)*u_scale, ',', v_tot(jx, jy, jz)*u_scale, ',',     &
            w_tot(jx, jy, jz)*u_scale, ',', theta_tot(jx, jy, jz)*T_scale 
    end do
    end do
    close(9)


    ! write(fname, '(a,i8.8,a)') '../output/vel_3D_check_', num, '.csv'
    ! open(9,file=fname)
    ! write(9, '(9a)') 'u', ',', 'v', ',', 'w', ',', 'theta', ',', 'PCon'
    ! do jz = 1, nzt-1
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     write(9, '(e15.7, a, e15.7, a, e15.7, a, e15.7, a, e15.7)')       &
    !         u_tot(jx, jy, jz), ',', v_tot(jx, jy, jz), ',',     &
    !         w_tot(jx, jy, jz), ',', theta_tot(jx, jy, jz), ',', &
    !         PCon_tot(jx, jy, jz, 1) 
    ! end do
    ! end do
    ! end do
    ! close(9)

    !  write(fname, '(a,i8.8,a)') '../output/vel_1D_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "z","u","v","w","p","theta"'
    !  write(9,*) 'zone i=', nzt-1, 'f=point'
    !  do jz = 1, nzt-1
    !      zz = (jz-1) * z_i / (nzt-1)
    !      write(9, '(4e15.7)') -zz, sum(u_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &      
    !                                sum(v_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &
    !                                sum(w_tot(:, :, jz))/float(nxt*nyt)*u_scale,         &
	! 			                   sum(p_tot(:, :, jz))/float(nxt*nyt),                 &
    !                                sum(theta_tot(:, :, jz))/float(nxt*nyt)*T_scale  
    !  end do
    !  close(9)
    
    !  write(fname, '(a,i8.8,a)') '../output/vel_xz_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "x","z","u","v","w","theta"'
    !  write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
    !  jy = nyt/2
    !  do jz = 1, nzt-1
    !  do jx = 1, nxt
    !      xx = (jx-1) * lx_tot * z_i / nxt
    !      zz = (jz-1) * z_i / (nzt-1)
    !      write(9, '(5e15.7)') xx, -zz, u_tot(jx, jy, jz)*u_scale,     &
    !                           v_tot(jx, jy, jz)*u_scale, w_tot(jx, jy, jz)*u_scale,                     &
    !                           theta_tot(jx, jy, jz)*T_scale 
    !  end do
    !  end do
    !  close(9)
    
    !  write(fname, '(a,i8.8,a)') '../output/vel_xy_check', num, '.dat'
    !  open(9,file=fname)
    !  write(9,*) 'variables = "x","y","u","v","w","theta"'
    !  write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
    !  jz = 2
    !  do jy = 1, nyt
    !  do jx = 1, nxt
    !      xx = (jx-1) * lx_tot * z_i / nxt
    !      yy = (jy-1) * ly_tot * z_i / nyt
    !      write(9, *) xx, yy, u_tot(jx, jy, jz)*u_scale,   &
    !                  v_tot(jx, jy, jz)*u_scale, w_tot(jx, jy, jz)*u_scale,                  &
    !                  theta_tot(jx, jy, jz)*T_scale 
    !  end do
    !  end do
    !  close(9)
    
    ! write(fname, '(a,i8.8,a)') '../output/vel_3D_check', num, '.dat'
    ! open(9,file=fname)
    ! write(9,*) 'variables = "x","y","z","u","v","w"'
    ! write(9,*) 'zone i=', nxt, 'j=', nyt, 'k=', nzt-1, 'f=point'
    ! do jz = 1, nzt-1
    ! do jy = 1, nyt
    ! do jx = 1, nxt
    !     xx = (jx-1) * lx_tot * z_i / nxt
    !     yy = (jy-1) * ly_tot * z_i / nyt
    !     zz = (jz-1) * z_i / (nzt-1)
    !     write(9, '(6e15.7)') xx, yy, zz, u_tot(jx, jy, jz), v_tot(jx, jy, jz), w_tot(jx, jy, jz) 
    ! end do
    ! end do
    ! end do
    ! close(9)
end if
! ---
end subroutine check_field
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine check_betascal(num)
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
use types, only: rprec
use param
use scalars_module, only: beta_scal
implicit none
! ---
integer, intent(in) :: num
real, dimension(ldx, nyt, nzt) :: beta_scal_tot

character(80) :: fname
integer :: jx, jy, jz
real(rprec) :: xx, yy, zz

! ---
call collect_data_3d(beta_scal, beta_scal_tot)

if ( rank == 0 ) then
     write(fname, '(a,i8.8,a)') '../output/beta_1D_check', num, '.dat'
     open(9,file=fname)
     write(9,*) 'variables = "z","beta"'
     write(9,*) 'zone i=', nzt-1, 'f=point'
     do jz = 1, nzt-1
         zz = (jz-1) * z_i / (nzt-1)
         write(9, '(4e15.4)') -zz, sum(beta_scal_tot(:, :, jz))/float(nxt*nyt)  
     end do
     close(9)
    
     write(fname, '(a,i8.8,a)') '../output/beta_xz_check', num, '.dat'
     open(9,file=fname)
     write(9,*) 'variables = "x","z","beta"'
     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
     jy = nyt/2
     do jz = 1, nzt-1
     do jx = 1, nxt
         xx = (jx-1) * lx_tot * z_i / nxt
         zz = (jz-1) * z_i / (nzt-1)
         write(9, '(5e15.4)') xx, -zz, beta_scal_tot(jx, jy, jz) 
     end do
     end do
     close(9)
    
     write(fname, '(a,i8.8,a)') '../output/beta_xy_check', num, '.dat'
     open(9,file=fname)
     write(9,*) 'variables = "x","y","beta"'
     write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
     jz = 4
     do jy = 1, nyt
     do jx = 1, nxt
         xx = (jx-1) * lx_tot * z_i / nxt
         yy = (jy-1) * ly_tot * z_i / nyt
         write(9, *) xx, yy, beta_scal_tot(jx, jy, jz)
     end do
     end do
     close(9)    
end if
! ---
end subroutine check_betascal
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
!subroutine check_pre(num)
!!-----------------------------------------------------------------------
!!   
!!-----------------------------------------------------------------------
!use types, only: rprec
!use param
!use sim_param, only: p
!implicit none
!! ---
!integer, intent(in) :: num
!real(rprec), dimension(ldx, nyt, 1:nzt) :: p_tot
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!
!! ---
!sendcounts = size (p(:,:,1:nz))
!recvcounts = size (p(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_gatherv(p(1,1,1), sendcounts, MPI_RPREC,           &
!                 p_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!    
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/pre_1D_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","p"'
!    write(9,*) 'zone i=', nzt-1, 'f=point'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') z, sum(p_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/pre_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","p"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(3e15.4)') x, z, p_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/pre_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","p"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, *) x, y, p_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_pre
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
!subroutine check_pg(num)
!!-----------------------------------------------------------------------
!!   
!!-----------------------------------------------------------------------
!use types, only: rprec
!use param
!use sim_param, only: dpdx, dpdy, dpdz
!implicit none
!! ---
!integer, intent(in) :: num
!real(rprec), dimension(ldx, nyt, 1:nzt) :: dpdx_tot
!real(rprec), dimension(ldx, nyt, 1:nzt) :: dpdy_tot
!real(rprec), dimension(ldx, nyt, 1:nzt) :: dpdz_tot
!
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!
!! ---
!sendcounts = size (dpdx(:,:,1:nz))
!recvcounts = size (dpdx(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_sendrecv(dpdx(1, 1, nz - 1), ldx*nyt, MPI_RPREC, up,   1, &
!                  dpdx(1, 1, 0),      ldx*nyt, MPI_RPREC, down, 1, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 2, &
!                  dpdx(1, 1, nz), ldx*nyt, MPI_RPREC, up,   2, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, nz - 1), ldx*nyt, MPI_RPREC, up,   3, &
!                  dpdx(1, 1, 0),      ldx*nyt, MPI_RPREC, down, 3, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 4, &
!                  dpdx(1, 1, nz), ldx*nyt, MPI_RPREC, up,   4, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, nz - 1), ldx*nyt, MPI_RPREC, up,   5, &
!                  dpdx(1, 1, 0),      ldx*nyt, MPI_RPREC, down, 5, &
!                  comm, status, ierr)
!call mpi_sendrecv(dpdx(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 6, &
!                  dpdx(1, 1, nz), ldx*nyt, MPI_RPREC, up,   6, &
!                  comm, status, ierr)
!
!call mpi_gatherv(dpdx(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdx_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(dpdy(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdy_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(dpdz(1,1,1), sendcounts, MPI_RPREC,           &
!                 dpdz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)                 
!    
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/pg_1D_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nzt-1, 'f=point'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') z, sum(dpdx_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(dpdy_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(dpdz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/pg_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(3e15.4)') x, z, dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz),  &
!                                   dpdz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/pg_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","dpdx","dpdy","dpdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = nz - 5
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, *) x, y, dpdx_tot(jx, jy, jz), dpdy_tot(jx, jy, jz),  &
!                          dpdz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_pg
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
! subroutine check_stokes_drift(num)
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
! use types, only: rprec
! use param
! use stokes_drift
! implicit none
! ! ---
! integer, intent(in) :: num
! real(rprec), dimension(0:nz) :: ust_abs
! real(rprec), dimension(0:nzt) :: ust_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs

! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: x, y, z
! ! ---

! ! ---
! end subroutine check_stokes_drift
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
! subroutine check_RHS(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num

! real, dimension(ldx, nyt, 1:nzt) :: RHSx_tot, RHSy_tot, RHSz_tot
! ! ---
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, yy, zz
! ! ---
! call collect_data_3d(RHSx, RHSx_tot)
! call collect_data_3d(RHSy, RHSy_tot)
! call collect_data_3d(RHSz, RHSz_tot)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/rhs_1D_check_', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "z","rhsx","rhsy","rhsz"'
!     do jz = 1, nzt-1
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(4e15.4)') -zz, sum(RHSx_tot(:, :, jz))/float(nxt*nyt),   & 
!                                 sum(RHSy_tot(:, :, jz))/float(nxt*nyt),   &
!                                 sum(RHSz_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)

!     write(fname, '(a,i8.8,a)') '../output/rhs_xz_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","rhsx","rhsy","rhsz"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(5e15.4)') xx, -zz, RHSx_tot(jx, jy, jz), RHSy_tot(jx, jy, jz),     &
!                                 RHSz_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)

!     write(fname, '(a,i8.8,a)') '../output/rhs_xy_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","y","rhsx","rhsy","rhsz"'
!     write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!     jz = 4
!     do jy = 1, nyt
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         yy = (jy-1) * ly_tot * z_i / nyt
!         write(9, *) xx, yy, RHSx_tot(jx, jy, jz), RHSy_tot(jx, jy, jz),     &
!                     RHSz_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_RHS
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
subroutine check_RHST(num)
!--------------------------------------------------------------------!
! 
!--------------------------------------------------------------------!
use types, only:rprec
use param
use scalars_module, only: RHS_T
implicit none
! ---
integer, intent(in) :: num

real, dimension(ldx, nyt, 1:nzt) :: RHST_tot
! ---
character(80) :: fname
integer :: jx, jy, jz
real(rprec) :: xx, yy, zz
! ---
call collect_data_3d(RHS_T, RHST_tot)

if ( rank == 0 ) then
    write(fname, '(a,i8.8,a)') '../output/rhsT_1D_check_', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "z","rhsT"'
    do jz = 1, nzt-1
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(4e15.4)') -zz, sum(RHST_tot(:, :, jz))/float(nxt*nyt)
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/rhsT_xz_check', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "x","z","rhsT"'
    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
    jy = nyt/2
    do jz = 1, nzt-1
    do jx = 1, nxt
        xx = (jx-1) * lx_tot * z_i / nxt
        zz = (jz-1) * z_i / (nzt-1)
        write(9, '(5e15.4)') xx, -zz, RHST_tot(jx, jy, jz)
    end do
    end do
    close(9)

    write(fname, '(a,i8.8,a)') '../output/rhsT_xy_check', num, '.dat'
    open(9,file=fname)
    write(9,*) 'variables = "x","y","rhsT"'
    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
    jz = 4
    do jy = 1, nyt
    do jx = 1, nxt
        xx = (jx-1) * lx_tot * z_i / nxt
        yy = (jy-1) * ly_tot * z_i / nyt
        write(9, *) xx, yy, RHST_tot(jx, jy, jz)
    end do
    end do
    close(9)
end if
! ---
end subroutine check_RHST
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
!subroutine check_divstress(num)
!!--------------------------------------------------------------------!
!! 
!!--------------------------------------------------------------------!
!use types, only:rprec
!use param
!use sim_param
!implicit none
!! ---
!integer, intent(in) :: num
!
!real(rprec), dimension(ldx, nyt, 1:nzt) :: divtx_tot, divty_tot, divtz_tot
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!! ---
!character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!! ---
!sendcounts = size (divtx(:,:,1:nz))
!recvcounts = size (divtx(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_gatherv(divtx(1,1,1), sendcounts, MPI_RPREC,           &
!                 divtx_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(divty(1,1,1), sendcounts, MPI_RPREC,           &
!                 divty_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(divtz(1,1,1), sendcounts, MPI_RPREC,           &
!                 divtz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/divstress_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","divtx","divty","divtz"'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(4e15.4)') z, sum(divtx_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(divty_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(divtz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/divstress_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","divtx","divty","divtz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(5e15.4)') x, z, divtx_tot(jx, jy, jz), divty_tot(jx, jy, jz),    &
!                                   divtz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/divstress_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(5e15.4)') x, y, divtx_tot(jx, jy, jz), divty_tot(jx, jy, jz),    &
!                                   divtz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_divstress
!!--------------------------------------------------------------------!
!!                                               
!!--------------------------------------------------------------------! 
!subroutine check_stress(num)
!!--------------------------------------------------------------------!
!! 
!!--------------------------------------------------------------------!
!use types, only:rprec
!use param
!use sim_param
!implicit none
!! ---
!integer, intent(in) :: num
!
!real(rprec), dimension(ldx, nyt, 1:nzt) :: txx_tot, txy_tot, txz_tot,     &
!                                             tyy_tot, tyz_tot, tzz_tot
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!! ---
!character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!! ---
!sendcounts = size (txx(:,:,1:nz))
!recvcounts = size (txx(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_sendrecv(txx(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 1,    &
!                  txx(1, 1, nz), ldx*nyt, MPI_RPREC, up,   1,    &
!                  comm, status, ierr)
!call mpi_sendrecv(txy(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 2,    &
!                  txy(1, 1, nz), ldx*nyt, MPI_RPREC, up,   2,    &
!                  comm, status, ierr)
!call mpi_sendrecv(tyy(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 3,    &
!                  tyy(1, 1, nz), ldx*nyt, MPI_RPREC, up,   3,    &
!                  comm, status, ierr)
!call mpi_sendrecv(tzz(1, 1, 1),  ldx*nyt, MPI_RPREC, down, 4,    &
!                  tzz(1, 1, nz), ldx*nyt, MPI_RPREC, up,   4,    &
!                  comm, status, ierr)
!
!call mpi_gatherv(txx(1,1,1), sendcounts, MPI_RPREC,           &
!                 txx_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(txy(1,1,1), sendcounts, MPI_RPREC,           &
!                 txy_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(txz(1,1,1), sendcounts, MPI_RPREC,           &
!                 txz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(tyy(1,1,1), sendcounts, MPI_RPREC,           &
!                 tyy_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(tyz(1,1,1), sendcounts, MPI_RPREC,           &
!                 tyz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_gatherv(tzz(1,1,1), sendcounts, MPI_RPREC,           &
!                 tzz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/stress_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","txx","txy","txz","tyy","tyz","tzz"'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(7e15.4)') z, sum(txx_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(txy_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(txz_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(tyy_tot(:, :, jz))/float(nxt*nyt),   & 
!                                sum(tyz_tot(:, :, jz))/float(nxt*nyt),   &
!                                sum(tzz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/stress_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(8e15.4)') x, z, txx_tot(jx, jy, jz), txy_tot(jx, jy, jz),    &
!                                   txz_tot(jx, jy, jz), tyy_tot(jx, jy, jz),    &
!                                   tyz_tot(jx, jy, jz), tzz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/stress_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","txx","txy","txz","tyy","tyz","tzz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(8e15.4)') x, y, txx_tot(jx, jy, jz), txy_tot(jx, jy, jz),    &
!                                   txz_tot(jx, jy, jz), tyy_tot(jx, jy, jz),    &
!                                   tyz_tot(jx, jy, jz), tzz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_stress
!!--------------------------------------------------------------------!
!!                                               
!!--------------------------------------------------------------------! 
!subroutine check_derivative(num)
!!--------------------------------------------------------------------!
!! 
!!--------------------------------------------------------------------!
!use types, only:rprec
!use param
!use sim_param
!implicit none
!! ---
!integer, intent(in) :: num
!
!real(rprec), dimension(ldx, nyt, 1:nzt) :: dwdz_tot
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!! ---
!character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
!character(80) :: fname
!integer :: jx, jy, jz
!real(rprec) :: x, y, z
!! ---
!sendcounts = size (dwdz(:,:,1:nz))
!recvcounts = size (dwdz(:,:,1:nz-1))
!displs = coord_of_rank * recvcounts
!
!call mpi_gatherv(dwdz(1,1,1), sendcounts, MPI_RPREC,           &
!                 dwdz_tot(1,1,1), sendcounts, displs,          &
!                 MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!if ( rank == 0 ) then
!    write(fname, '(a,i8.8,a)') '../output/dwdz_1D_check_', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "z","dwdz"'
!    do jz = 1, nzt-1
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(2e15.4)') z, sum(dwdz_tot(:, :, jz))/float(nxt*nyt)
!    end do
!    close(9)
!
!    write(fname, '(a,i8.8,a)') '../output/dwdz_xz_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","z","dwdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!    jy = nyt/2
!    do jz = 1, nzt-1
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        z = (jz-1) * z_i / (nzt-1)
!        write(9, '(3e15.4)') x, z, dwdz_tot(jx, jy, jz)
!    end do
!    end do
!    close(9)
!    
!    write(fname, '(a,i8.8,a)') '../output/dwdz_xy_check', num, '.dat'
!    open(9,file=fname)
!    write(9,*) 'variables = "x","y","dwdz"'
!    write(9,*) 'zone i=', nxt, 'j=', nyt, 'f=point'
!    jz = 4
!    do jy = 1, nyt
!    do jx = 1, nxt
!        x = (jx-1) * lx_tot * z_i / nxt
!        y = (jy-1) * ly_tot * z_i / nyt
!        write(9, '(8e15.4)') x, y, dwdz_tot(jx, jy, jz) 
!    end do
!    end do
!    close(9)
!end if
!! ---
!end subroutine check_derivative
!!--------------------------------------------------------------------!
!!                                               
!!--------------------------------------------------------------------!
! subroutine check_divt(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! ! ---
! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/divt_check_', num, '.csv'
!     open(9,file=fname)
!     write(9, '(5a)') 'divtx', ',', 'divty', ',', 'divtxz'
!     do jz = 1, nz-1
!     do jy = 1, nyt
!     do jx = 1, nxt
!         write(9, fmt) divtx(jx, jy, jz), ',', divty(jx, jy, jz), ',', divtz(jx, jy, jz)
!     end do
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_divt
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------! 
! subroutine check_coriolis(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use types, only:rprec
! use param
! use sim_param
! implicit none
! ! ---
! integer, intent(in) :: num
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: xx, zz
! ! ---
! if ( rank == 0 ) then
!     ! write(fname, '(a,i8.8,a)') '../output/coriolis_check_', num, '.csv'
!     ! open(9,file=fname)
!     ! write(9, '(5a)') 'coriol_v', ',', 'coriol_u'
!     ! do jz = 1, nz-1
!     ! do jy = 1, nyt
!     ! do jx = 1, nxt
!     !     write(9, fmt) coriol*v(jx, jy, jz), ',', coriol*u(jx, jy, jz)
!     ! end do
!     ! end do
!     ! end do
!     ! close(9)

!     write(fname, '(a,i8.8,a)') '../output/coriolis_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","fv","fu"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         xx = (jx-1) * lx_tot * z_i / nxt
!         zz = (jz-1) * z_i / (nzt-1)
!         write(9, '(4e15.4)') xx, zz, coriol*v(jx, jy, jz), coriol*u(jx, jy, jz) 
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_coriolis
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------! 
! subroutine check_Nut(num)
! !--------------------------------------------------------------------!
! ! 
! !--------------------------------------------------------------------!
! use param
! use sgsmodule, only: Cs_opt2, Nu_t
! implicit none
! ! ---
! integer, intent(in) :: num
! real(rprec), dimension(ldx, nyt, 1:nzt) :: Nu_t_tot
! integer, dimension(npz) :: sendcounts, recvcounts, displs
! ! ---
! character(*), parameter :: fmt = '(e15.4, a, e15.4, a, e15.4, a, e15.4)'
! character(80) :: fname
! integer :: jx, jy, jz
! real(rprec) :: x, y, z
! ! ---
! sendcounts = size (Nu_t(:,:,1:nz))
! recvcounts = size (Nu_t(:,:,1:nz-1))
! displs = coord_of_rank * recvcounts

! call mpi_gatherv(Nu_t(1,1,1), sendcounts, MPI_RPREC,           &
!                  Nu_t_tot(1,1,1), sendcounts, displs,          &
!                  MPI_RPREC, rank_of_coord(0), comm, ierr)

! if ( rank == 0 ) then
!     write(fname, '(a,i8.8,a)') '../output/Nut_1D_check_', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "z","Nu_t"'
!     write(9,*) 'zone i=', nzt-1, 'f=point'
!     do jz = 1, nzt-1
!         z = (jz-1) * z_i / (nzt-1)
!         write(9, '(2e15.4)') z, sum(Nu_t_tot(:, :, jz))/float(nxt*nyt)
!     end do
!     close(9)
    
!     write(fname, '(a,i8.8,a)') '../output/Nut_xz_check', num, '.dat'
!     open(9,file=fname)
!     write(9,*) 'variables = "x","z","Nu_t"'
!     write(9,*) 'zone i=', nxt, 'j=', nzt-1, 'f=point'
!     jy = nyt/2
!     do jz = 1, nzt-1
!     do jx = 1, nxt
!         x = (jx-1) * lx_tot * z_i / nxt
!         z = (jz-1) * z_i / (nzt-1)
!         write(9, '(3e15.4)') x, z, Nu_t_tot(jx, jy, jz)
!     end do
!     end do
!     close(9)
! end if
! ! ---
! end subroutine check_Nut
!--------------------------------------------------------------------!
                                               
!--------------------------------------------------------------------! 