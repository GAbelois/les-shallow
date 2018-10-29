!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(in)  :: a_local
real, dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))
call mpi_gather(a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data2_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nynpy, nz), intent(in)  :: a_local
real, dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, nz, 0:np-1))
call mpi_gather(a_local, ldx*nynpy*nz, mpi_rprec,   &
                temp,    ldx*nynpy*nz, mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do

end do

deallocate(temp)
! ---
end subroutine collect_data2_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data3_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nxt, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(nxt, nynpy, nz), intent(in)  :: a_local
real, dimension(nxt, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(nxt, nynpy, nz, 0:np-1))
call mpi_gather(a_local, nxt*nynpy*nz, mpi_rprec,   &
                temp,    nxt*nynpy*nz, mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, nxt
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do

end do

deallocate(temp)
! ---
end subroutine collect_data3_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data4_3d(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                    idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(in)  :: a_local
real, dimension(ldx, nyt, nzt),  intent(out) :: a_global
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, 0:nynpy+1, 0:nz, 0:np-1))
call mpi_gather(a_local, ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                temp,    ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                0, comm, ierr)

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz) = temp(jx, jy, jz, ip) 
    end do
    end do
    end do
end do

deallocate(temp)
! ---
end subroutine collect_data4_3d
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_z(a_local, a_global)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: nz, nzt, npz, idz_col, mpi_rprec, comm_col, ierr
implicit none
! ---
real(kind=rprec), dimension(nz-1),  intent(in)  :: a_local
real(kind=rprec), dimension(nzt-1), intent(out) :: a_global
! ---
real(kind=rprec), dimension(nz-1, 0:npz-1) :: temp
integer :: ip, ipz, jz
! ---
call mpi_allgather(a_local, nz-1, mpi_rprec, temp, nz-1, mpi_rprec, comm_col, ierr)

do ipz = 0, npz - 1
    ip = idz_col(ipz)
    do jz = 1, nz - 1
        a_global(ip*(nz-1)+jz) = temp(jz, ipz)
    end do
end do
! ---
end subroutine collect_data_z
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                 a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data2(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                    idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real(kind=rprec), dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*nz, mpi_rprec,   &
                 a_local, ldx*nynpy*nz, mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine scatter_data3(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param, only: ldx, nynpy, nyt, nz, nzt, np,    &
                 idy, idz, mpi_rprec, comm, ierr
implicit none
! ---
real, dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, nynpy, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, nynpy, 0:nz, 0:np-1))

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*nynpy*(nz+1), mpi_rprec,   &
                 a_local, ldx*nynpy*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

deallocate(temp)
! ---
end subroutine scatter_data3
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine scatter_data4(a_global, a_local)
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
use types, only: rprec
use param
implicit none
! ---
real, dimension(ldx, nyt, nzt),  intent(in)  :: a_global
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:, :, :, :), allocatable :: temp
real(kind=rprec), dimension(:,:), allocatable :: temp_xz0, temp_xz1, temp_xz2, temp_xz3
integer :: ip, ipy, ipz, jx, jy, jz
! ---
allocate(temp(ldx, 0:nynpy+1, 0:nz, 0:np-1))
allocate(temp_xz0(ldx, 0:nz), temp_xz1(ldx, 0:nz), temp_xz2(ldx, 0:nz), temp_xz3(ldx, 0:nz))
temp = 0._rprec
temp_xz0 = 0._rprec
temp_xz1 = 0._rprec
temp_xz2 = 0._rprec
temp_xz3 = 0._rprec

do ip = 0, np - 1
    ipy = idy(ip)
    ipz = idz(ip)
    do jz = 1, nz
    do jy = 1, nynpy
    do jx = 1, ldx
        temp(jx, jy, jz, ip) = a_global(jx,ipy*nynpy+jy,ipz*(nz-1)+jz)
    end do
    end do
    end do
end do

call mpi_scatter(temp,    ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                 a_local, ldx*(nynpy+2)*(nz+1), mpi_rprec,   &
                 0, comm, ierr)

temp_xz0(:, :) = a_local(:, 1, :)
call mpi_sendrecv(temp_xz0(1, 1), ldx*(nz+1), MPI_RPREC, south, 1,       &
                  temp_xz1(1, 1), ldx*(nz+1), MPI_RPREC, north, 1,       &
                  comm, status, ierr)
a_local(:,nynpy+1,:) = temp_xz1(:, :)


temp_xz2(:, :) = a_local(:, nynpy, :)
call mpi_sendrecv(temp_xz2(1, 1), ldx*(nz+1), MPI_RPREC, south, 2,       &
                  temp_xz3(1, 1), ldx*(nz+1), MPI_RPREC, north, 2,       &
                  comm, status, ierr)
a_local(:,0,:) = temp_xz3(:, :)                  
                  
deallocate(temp, temp_xz0, temp_xz1, temp_xz2, temp_xz3)
! ---
end subroutine scatter_data4
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine collect_data_xh(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nxhnpy, lhx, nyt, npy, coordz, idy_ver,  &
                 comm_ver, mpi_cprec, ierr
implicit none
! ---
complex(rprec), dimension(nxhnpy, nyt), intent(in) :: in
complex(rprec), dimension(lhx, nyt),   intent(out) :: out
! ---
complex(rprec), dimension(nxhnpy, nyt, 0:npy-1) :: tmp
integer :: ip, ipy, jy, num
! ---
num = nxhnpy * nyt
call mpi_allgather(in, num, mpi_cprec, tmp, num, mpi_cprec, comm_ver, ierr)

do ipy = 0, npy-1
    ip = idy_ver(ipy)
    do jy = 1, nxhnpy
        out(ip*(nxhnpy-1)+jy, :) = tmp(jy, :, ipy)
    end do
end do
! ---
end subroutine collect_data_xh
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine collect_data_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, npy, coordz, idy_ver, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nynpy), intent(in) :: in
real(rprec), dimension(nyt),  intent(out) :: out
! ---
real(rprec), dimension(nynpy, 0:npy-1) :: tmp
integer :: ip, ipy, jy
! ---
call mpi_allgather(in, nynpy, mpi_rprec, tmp, nynpy, mpi_rprec, comm_ver, ierr)

do ipy = 0, npy-1
    ip = idy_ver(ipy)
    do jy = 1, nynpy
        out(ip*nynpy+jy) = tmp(jy, ipy)
    end do
end do
! ---
end subroutine collect_data_y 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine scatter_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, npy, idy_ver, comm_ver, mpi_rprec, ierr
implicit none
! ---
real(rprec), dimension(nyt),   intent(in)  :: in
real(rprec), dimension(nynpy), intent(out) :: out
! ---
real(rprec), dimension(nynpy, 0:npy-1) :: tmp
integer :: ipy, j, ind
! ---
do ipy = 0, npy-1
    do j = 1, nynpy
        !ind = ipy*nynpy + j
        !tmp(j,rank_ver(ipy)) = in(ind)
        ind = idy_ver(ipy)*nynpy + j
        tmp(j,ipy) = in(ind)    
    end do
end do

call mpi_scatter(tmp, nynpy, mpi_rprec,     &
                 out, nynpy, mpi_rprec, 0, comm_ver, ierr)
! ---
end subroutine scatter_y
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine data_exchange()
!-----------------------------------------------------------------------
!   exchange data between neighboring zones
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
integer :: ipcon, jx, jy, jz
real :: temp_w
! --- synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz
call mpi_sendrecv(u(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   1,  &
                  u(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 1,  &
                  comm, status, ierr)
call mpi_sendrecv(v(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   2,  &
                  v(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 2,  &
                  comm, status, ierr)
call mpi_sendrecv(w(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   3,  &
                  w(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 3,  &
                  comm, status, ierr)
call mpi_sendrecv(u(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 4,  &
                  u(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   4,  &
                  comm, status, ierr)
call mpi_sendrecv(v(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 5,  &
                  v(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   5,  &
                  comm, status, ierr)
call mpi_sendrecv(w(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 6,  &
                  w(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   6,  &
                  comm, status, ierr)
if (theta_flag) then
    call mpi_sendrecv(theta(1, 1, nz - 1), ldx*nynpy, MPI_RPREC, up,   7,  &
                      theta(1, 1, 0),      ldx*nynpy, MPI_RPREC, down, 7,  &
                      comm, status, ierr)
    call mpi_sendrecv(theta(1, 1, 1),  ldx*nynpy, MPI_RPREC, down, 8,  &
                      theta(1, 1, nz), ldx*nynpy, MPI_RPREC, up,   8,  &
                      comm, status, ierr)
end if

if (PCon_FLAG) then
    do ipcon = 1, npcon
        call mpi_sendrecv(PCon(1, 0, nz - 1, ipcon), ldx*(nynpy+1), MPI_RPREC, up,   1,    &
                          PCon(1, 0, 0, ipcon),      ldx*(nynpy+1), MPI_RPREC, down, 1,    &
                          comm, status, ierr)
    end do
end if

if ( coordz == 0 ) then
    !--set 0-level velocities to BOGUS
    u(:, :, 0) = BOGUS
    v(:, :, 0) = BOGUS
    w(:, :, 0) = BOGUS
    if (theta_flag) theta(:, :, 0) = BOGUS
else if ( coordz == npz - 1 ) then
    u(:, :, nz) = BOGUS
    v(:, :, nz) = BOGUS
    !w(:, :, nz) = BOGUS
    if (theta_flag) theta(:, :, nz) = BOGUS
end if
! ---
end subroutine data_exchange
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine data_exchange2(a_local)
!-----------------------------------------------------------------------
!   exchange data for PCon
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, 0:nz), intent(out) :: a_local
! ---
real(kind=rprec), dimension(:,:), allocatable :: temp_xz0, temp_xz1, temp_xz2, temp_xz3

integer :: ipcon, jx, jy, jz
! ---
allocate(temp_xz0(ldx, 0:nz), temp_xz1(ldx, 0:nz), temp_xz2(ldx, 0:nz), temp_xz3(ldx, 0:nz))
temp_xz0 = 0._rprec
temp_xz1 = 0._rprec
temp_xz2 = 0._rprec
temp_xz3 = 0._rprec

! ---
temp_xz0(:, :) = a_local(:, 1, :)
call mpi_sendrecv(temp_xz0(1, 0), ldx*(nz+1), MPI_RPREC, south, 1,       &
                  temp_xz1(1, 0), ldx*(nz+1), MPI_RPREC, north, 1,       &
                  comm, status, ierr)
a_local(:,nynpy+1,:) = temp_xz1(:, :)


temp_xz2(:, :) = a_local(:, nynpy, :)
call mpi_sendrecv(temp_xz2(1, 0), ldx*(nz+1), MPI_RPREC, north, 2,       &
                  temp_xz3(1, 0), ldx*(nz+1), MPI_RPREC, south, 2,       &
                  comm, status, ierr)
a_local(:,0,:) = temp_xz3(:, :) 

! --- synchronize the overlapping parts 0 <-> nz-1 and 1 <-> nz
call mpi_sendrecv(a_local(1, 0, nz - 1), ldx*(nynpy+2), MPI_RPREC, up,   3,    &
                  a_local(1, 0, 0),      ldx*(nynpy+2), MPI_RPREC, down, 3,    &
                  comm, status, ierr)
call mpi_sendrecv(a_local(1, 0, 1),      ldx*(nynpy+2), MPI_RPREC, down, 4,    &
                  a_local(1, 0, nz),     ldx*(nynpy+2), MPI_RPREC, up,   4,    &
                  comm, status, ierr)
! ---
deallocate(temp_xz0, temp_xz1, temp_xz2, temp_xz3)
! ---
end subroutine data_exchange2
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine data_exchange_y(a_local, y_per)
!-----------------------------------------------------------------------
!   exchange data for PCon
!-----------------------------------------------------------------------
use param
use sim_param
implicit none
! ---
real(kind=rprec), dimension(ldx, 0:nynpy+1, nz), intent(out) :: a_local
logical, intent(in) :: y_per 
! ---
real(kind=rprec), dimension(ldx,nz) :: temp0, temp1, temp0m, temp1m
integer :: sou, nor
! ---
temp0  = 0._rprec
temp1  = 0._rprec
temp0m = 0._rprec
temp1m = 0._rprec

! ---
temp1(:, :)  = a_local(:, 1, :)
temp1m(:, :) = a_local(:, nynpy, :)

if ( idy(rank) > 0 ) then
    sou = rank - npz
else 
    sou = mpi_proc_null
end if

if ( idy(rank) < npy-1 ) then
    nor = rank + npz
else
    nor = mpi_proc_null
end if

call mpi_sendrecv(temp1(:,:),  ldx*nz, MPI_RPREC, sou,  1,  &
                  temp0m(:,:), ldx*nz, MPI_RPREC, nor,  1,  &
                  comm, status, ierr)
call mpi_sendrecv(temp1m(:,:), ldx*nz, MPI_RPREC, nor,  2,  &
                  temp0(:,:),  ldx*nz, MPI_RPREC, sou,  2,  &
                  comm, status, ierr)
! --- periodic condition
if ( y_per ) then
    if ( idy(rank) == 0 ) then
        sou = rank + (npy-1)*npz
    else
        sou = mpi_proc_null
    end if

    if ( idy(rank) == npy-1 ) then
        nor = rank - (npy-1)*npz
    else
        nor = mpi_proc_null
    end if

    call mpi_sendrecv(temp1(:,:),  ldx*nz, MPI_RPREC, sou,  1,  &
                      temp0m(:,:), ldx*nz, MPI_RPREC, nor,  1,  &
                      comm, status, ierr)
    call mpi_sendrecv(temp1m(:,:), ldx*nz, MPI_RPREC, nor,  2,  &
                      temp0(:,:),  ldx*nz, MPI_RPREC, sou,  2,  &
                      comm, status, ierr)
end if

if ( idy(rank) > 0 .or. y_per ) then
    a_local(:,0,:) = temp0(:, :)
end if

if ( idy(rank) < npy-1 .or. y_per ) then
    a_local(:,nynpy+1,:) = temp0m(:, :)
end if

! ---
end subroutine data_exchange_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine avg_y(in, out)
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
use types, only: rprec
use param, only: nynpy, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nynpy), intent(in)  :: in
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
integer :: ipy, j, ind
! ---
tmp = 0._rprec
do j = 1, nynpy
    tmp = tmp + in(j)
end do
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)

out = out / nyt
! ---
end subroutine avg_y
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
subroutine calc_hor_avg(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(ldx, nynpy, 0:nz), intent(in)  :: in
real(rprec), dimension(0:nz-1), intent(out) :: out
! ---
real(rprec), dimension(0:nz-1) :: tmp
integer :: jz
! ---
do jz = 0, nz-1
    tmp(jz) = sum(in(1:nxt, 1:nynpy, jz)) / (nxt * nyt)
end do
call mpi_allreduce(tmp, out, nz, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
subroutine calc_hor_avg2(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(ldx, nynpy, nz), intent(in)  :: in
real(rprec), dimension(nz-1), intent(out) :: out
! ---
real(rprec), dimension(nz-1) :: tmp
integer :: jz
! ---
do jz = 1, nz-1
    tmp(jz) = sum (in(1:nxt, 1:nynpy, jz)) / (nxt * nyt)
end do
call mpi_allreduce(tmp, out, nz-1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg2
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine calc_hor_avg3(in, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: ldx, nynpy, nz, nxt, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nxt, nynpy), intent(in)  :: in
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
integer :: jz
! ---
tmp = sum(in(1:nxt, 1:nynpy)) / (nxt * nyt)
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_hor_avg3
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine calc_rms(in1, in2, out)    
!-----------------------------------------------------------------------
!   calculate the horizontal-average value of the 3d field 
!-----------------------------------------------------------------------     
use types, only: rprec
use param, only: nxt, nynpy, nyt, comm_ver, mpi_rprec, ierr
use mpi
implicit none
! ---
real(rprec), dimension(nxt, nynpy), intent(in)  :: in1, in2
real(rprec), intent(out) :: out
! ---
real(rprec) :: tmp
! ---
tmp = sum (in1*in2) / (nxt * nyt)
call mpi_allreduce(tmp, out, 1, mpi_rprec, mpi_sum, comm_ver, ierr)
! ---    
end subroutine calc_rms
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
!subroutine collocate_MPI_averages_N(avg_var_proc, avg_var_tot_domain, file_ind, filename_str)
!!-----------------------------------------------------------------------
!!--The following subroutine does the collocation of the MPI arrays for
!! averaging in avgslice and scalar_slice (in scalars_module2.f90)
!!-----------------------------------------------------------------------
!use types, only: rprec
!use param
!use intermediate, only: dim1_size, dim2_size, dim1_global, dim2_global
!implicit none
!! ---
!real(kind=rprec), dimension(dim1_size, dim2_size)     :: avg_var_proc
!real(kind=rprec), dimension(dim1_global, dim2_global) :: avg_var_tot_domain
!integer, intent(in) :: file_ind
!character(*), intent(in) :: filename_str
!    
!integer :: ind1, ind2, jx
!character(len=256) :: local_filename
!character(*), parameter :: fmt_5168 = "(1400(E14.5))"
!
!local_filename = path//'result/aver_'//trim(filename_str)//'.out'
!avg_var_tot_domain = 0._rprec
!   
!if ( average_dim_num .eq. 1 ) then
!    do jx = 1, dim1_size
!        call collect_data_z(avg_var_proc(jx, 1:nz-1), avg_var_tot_domain(jx, 1:nzt-1))    
!    end do
!else if ( average_dim_num .eq. 2 ) then 
!    call collect_data_z(avg_var_proc(1:nz-1,1), avg_var_tot_domain(1:nzt-1,1))
!end if
!
!if ( rank == 0 ) then
!    open (file_ind, file=trim(local_filename), status="unknown", position="append")
!
!    if (average_dim_num .eq. 1) then
!        do ind2 = 1, nzt - 1
!        do ind1 = 1, nxt
!            if (abs(avg_var_tot_domain(ind1, ind2)) .lt. TINYS) then
!                avg_var_tot_domain(ind1, ind2) = 0._rprec
!            end if
!        end do
!        end do
!      
!        do ind2 = 1, nzt - 1
!            write (file_ind, fmt_5168) jt*dt, (avg_var_tot_domain(ind1, ind2), ind1=1, nxt)
!        end do
!    else if (average_dim_num .eq. 2) then
!        do ind1 = 1, nzt - 1
!            if (abs(avg_var_tot_domain(ind1, 1)) .lt. TINYS) then
!                avg_var_tot_domain(ind1, 1) = 0._rprec
!            end if
!        end do
!        write (file_ind, fmt_5168) jt*dt, (avg_var_tot_domain(ind1, 1), ind1=1, nzt - 1)
!    end if
!    
!    close (file_ind)
!end if
!! ---
!end subroutine collocate_MPI_averages_N
!!-----------------------------------------------------------------------
!!   
!!----------------------------------------------------------------------- 
