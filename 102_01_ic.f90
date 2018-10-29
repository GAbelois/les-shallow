!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
subroutine ic()
!--------------------------------------------------------------------!
! Log profile that is modified to flatten at z=z_i
!--------------------------------------------------------------------!
use types, only:rprec
use param
use sim_param, only:u, v, w
use bottombc
implicit none
! ---
real(kind=rprec), parameter :: alpha = 0.5_rprec
integer :: jx, jy, jz, seed, jz_abs

real(kind=rprec), dimension(nz) :: ubar, vbar
real(kind=rprec) :: rms, noise, arg, arg2
real(kind=rprec) :: z, w_star, ran3
! ---
if ((inflow) .and. (.not. read_inflow_file)) then !--no turbulence
    u = face_avg
    v = 0._rprec
    w = 0._rprec
else
    !BC uncorelate w_star from wt_s by Bicheng Chen
    w_star = (9.81_rprec/T_init*0.06_rprec*z_i)**(1._rprec/3._rprec)
    !w_star=(9.81_rprec/T_init*wt_s*z_i)**(1._rprec/3._rprec)
    !      T_star=wt_s/w_star
    !      q_star=T_star
    !write(*,*) 'Modified Log Profile for IC'
    
    if ( ocean_flag ) then
        call ekman_layer(ubar, vbar)
        ! do jz = 1, nz-1
        !     if ( ubc_mom == 'wall' ) then               ! shallow ocean flow
        !         select case (walltype)
        !         case('smooth')
        !             arg2 = (lz_tot - zuvp(jz)) / (nu_molec/u_scale/z_i)
        !             arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
        !         case('rough')
        !             arg2 = (lz_tot - zuvp(jz)) /(zo1/z_i)
        !             arg = (1._rprec/vonk)*log(arg2)
        !             !arg = 2._rprec*((lz_tot - zuvp(jz))/lz_tot)**0.3
        !         case default
        !             write(*,*) 'invalid wall type'
        !             stop
        !         end select
        !     else if ( ubc_mom == 'stress free' ) then   ! deep ocean flow
        !         ! waits to be implemented
        !         arg2 = zuvp(jz) /(zo1/z_i)
        !         arg = 1._rprec
        !         !arg = (1._rprec/vonk)*log(arg2)
        !     end if
        
        !     if ( coriolis_forcing ) then
        !         ubar(jz) = ug
        !         vbar(jz) = vg
        !     else
        !         ubar(jz) = arg
        !         vbar(jz) = 0._rprec
        !     end if
        
        !     !if ((coriolis_forcing) .and. (zuvp(jz) .gt. (.5_rprec*z_i))) ubar(jz) = ug 
        ! end do
    else   
        do jz = 1, nz-1 
            if ( ubc_mom == 'wall' ) then   ! channel flow with bottom and top walls
                select case (walltype)      ! bottom boundary
                case('smooth')
                    if ( zuvp(jz) .le. lz_tot/2 ) then
                        arg2 = zuvp(jz) / (nu_molec/u_scale/z_i)
                    else
                        arg2 = (lz_tot-zuvp(jz)) / (nu_molec/u_scale/z_i)
                    end if
                    arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
                case('rough')
                    ! IC in equilibrium with rough surface (rough dominates in effective zo)
                    ! arg2=z/(sum(zo)/float(nxt*nyt))
                    if ( zuvp(jz) .le. lz_tot/2 ) then
                        arg2 = zuvp(jz)/(zo1/z_i)
                    else
                        arg2 = (lz_tot-zuvp(jz))/(zo1/z_i)
                    end if
                    arg = (1._rprec/vonk)*log(arg2)
                case default
                    write(*,*) 'invalid wall type'
                    stop
                end select
            else if ( ubc_mom == 'stress free' ) then   ! atmospheric boundary layer flow
                select case (walltype)      ! bottom boundary
                case('smooth')
                    arg2 = zuvp(jz) / (nu_molec/u_scale/z_i)
                    arg = (1._rprec/vonk)*log(arg2) + B !-1./(2.*vonk*z_i*z_i)*z*z
                case('rough')
                    ! IC in equilibrium with rough surface (rough dominates in effective zo)
                    ! arg2=z/(sum(zo)/float(nxt*nyt))
                    arg2 = zuvp(jz)/(zo1/z_i)
                    arg = (1._rprec/vonk)*log(arg2)
                case default
                    write(*,*) 'invalid wall type'
                    stop
                end select
            end if
                
            if ( coriolis_forcing ) then
                ubar(jz) = ug
                vbar(jz) = vg
                ! Note that ug and vg have already been non-dimensionalized in param.f90
                !           ubar(jz)=arg/30._rprec
            else
                ubar(jz) = arg
                vbar(jz) = 0._rprec
            end if

            if ((coriolis_forcing) .and. (zuvp(jz) .gt. (.5_rprec*z_i))) ubar(jz) = ug 
        end do
    end if
    
    rms = 1._rprec
    do jz = 1, nz-1
        jz_abs = coordz*(nz - 1) + jz
        seed = -80 - jz_abs !--trying to make consistent init for MPI
        do jy = 1, nynpy
        do jx = 1, nxt
            !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
            !...Taking std dev of vel as 1 at all heights
            if (zuvp(jz) .le. lz_tot) then
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                u(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz)) + ubar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                !v(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz)) + vbar(jz)
                v(jx, jy, jz) = vbar(jz)
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                !w(jx, jy, jz) = alpha*noise*abs(0.5_rprec*lz_tot - zuvp(jz))
                w(jx, jy, jz) = 0._rprec
            else
                noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
                u(jx, jy, jz) = noise*w_star/u_scale*.01_rprec + ubar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                v(jx, jy, jz) = noise*w_star/u_scale*.01_rprec + vbar(jz)
                noise = rms/.289_rprec*(ran3(seed) - 0.5_rprec)
                w(jx, jy, jz) = noise*w_star/u_scale*.01_rprec
            end if
        end do
        end do
    end do

    !...BC for W
    if ( coordz == 0 ) then
        w(1:nxt, 1:nynpy, 1) = 0._rprec
    end if
    if ( coordz == npz - 1 ) then
        w(1:nxt, 1:nynpy, nz) = 0._rprec
    endif

    !...BC for U, V
    if ( coordz == npz - 1 ) then
        if ( ubc_mom == 'wall' ) then
            u(1:nxt, 1:nynpy, nz) = BOGUS
            v(1:nxt, 1:nynpy, nz) = BOGUS
        else if ( ubc_mom == 'stress free' ) then
            u(1:nxt, 1:nynpy, nz) = u(1:nxt, 1:nynpy, nz - 1)
            v(1:nxt, 1:nynpy, nz) = v(1:nxt, 1:nynpy, nz - 1)
        end if
    end if
end if
! ---
end subroutine ic
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
subroutine ic_dns()
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!    
use types, only:rprec
use param
use sim_param, only:u, v, w
implicit none
! ---  
real(kind=rprec), dimension(nz)::ubar
real(kind=rprec) :: rms, temp
integer :: jx, jy, jz, seed
real(kind=rprec) :: z

real(rprec), external :: ran3
! ---
if (inflow) then
    ! uniform flow case:
    ubar = face_avg
else
    seed = -112

    ! calculate height of first uvp point in wall units
    ! lets do a laminar case (?)
    ! *** cyan we shouldn't impose laminar flow profile using u_scale as velocity scale
    ! ***      because u_scale in laminar and turbulent scenarios are completely different
    do jz = 1, nz
        ubar(jz) = (u_scale*z_i/nu_molec)*zuvp(jz)*(1._rprec - .5_rprec*zuvp(jz)) ! non-dimensional
    end do
end if

! rms=0.0001 seems to work in some cases
! the "default" rms of a unif variable is 0.289
rms = 0.2_rprec
do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    u(jx, jy, jz) = ubar(jz) + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
    v(jx, jy, jz) = 0._rprec + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
    w(jx, jy, jz) = 0._rprec + (rms/.289_rprec)*(ran3(seed) - .5_rprec)/u_scale
end do
end do
end do

! make sure w-mean is 0
temp = 0._rprec
do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    temp = temp + w(jx, jy, jz)
end do
end do
end do
temp = temp/(nxt*nynpy*nz)

do jz = 1, nz
do jy = 1, nynpy
do jx = 1, nxt
    w(jx, jy, jz) = w(jx, jy, jz) - temp
end do
end do
end do

w(:, :, 1)  = 0._rprec
w(:, :, nz) = 0._rprec
u(:, :, nz) = u(:, :, nz - 1)
v(:, :, nz) = v(:, :, nz - 1)
! ---
end subroutine ic_dns
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine ic_scal()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param, only:u, v, w, theta, q
use bottombc
implicit none
! ---

real(kind=rprec), dimension(nz) :: ubar, vbar, wbar
real(kind=rprec) :: rms, noise, arg, arg2, theta_mean, ustar_ratio
real(kind=rprec) :: z, w_star, T_star, q_star, ran3, bv_freq
real(kind=rprec) :: z_turb_limit, perturb_height_factor, z_inv
integer :: jx, jy, jz, seed, jz_abs

character(*), parameter :: fmt_7781 = "('z, ubar, vbar, wbar, Tbar:',5(1x, F9.4))"
! ---
if ( lbc .eq. 0 ) then
    theta_mean = T_s_min - 0.001*T_s_min !T_s_min is dimensional while T_s is non-dimensional
else 
    theta_mean = T_init
end if

if (wt_s .lt. 0._rprec) then
    perturb_height_factor = 0.30_rprec
    z_inv = 0.30_rprec*z_i
else
    perturb_height_factor = 0.3_rprec
    z_inv = 0.57_rprec*z_i

    if (ocean_flag) then
        z_inv = prop_mixed*z_i
    end if
end if
z_turb_limit = perturb_height_factor*z_i

if (wt_s .eq. 0.0_rprec) then
    ! w_star is of O(1) with z_i=500 and wt_s=0.06 for atmos
    ! For ocean, ~5e-5 is better
    w_star = (g/theta_mean*5.3691d-5*z_i)**(1._rprec/3._rprec)
    T_star = 5.3691d-5/w_star
    q_star = T_star
else
    w_star = sign((g/theta_mean*abs(wt_s)*z_i)**(1._rprec/3._rprec), wt_s)
    T_star = wt_s/w_star
    q_star = T_star
end if

if (ocean_flag) then
    !call ekman_layer(ubar, vbar)
    if ( ubc_mom == 'wall' ) then
        call init_profile(ubar, vbar, z_inv)
        ! do jz = 1, nz-1
        !     !if ( zuvp(jz) .ge. lz_tot - ug_dim/bv_freq/z_i ) then
        !     ! select case (walltype)
        !     ! case ('smooth') 
        !     !     arg2 = (lz_tot-zuvp(jz)) / (nu_molec/u_scale/z_i)
        !     !     arg = (1._rprec/vonk)*log(arg2) + B
        !     ! case ('rough')
        !     !     arg2 = (lz_tot-zuvp(jz))/(zo1/z_i)
        !     !     arg = (1._rprec/vonk)*log(arg2)
        !     ! case default
        !     !     write(*,*) 'invalid wall type'
        !     !     stop
        !     ! end select
        !     ! !end if
        !     ! ubar(jz) = arg
        !     ! vbar(jz) = 0._rprec
        !     ! if ( zuvp(jz) .le. lz_tot/2._rprec ) then
        !     !     ubar(jz) = ug*(1+cosh(16._rprec*pi*(-zuvp(jz)+lz_tot))/     &
        !     !                      2._rprec/(sinh(8._rprec*pi*lz_tot))**2)
        !     ! else
        !     !     ubar(jz) = ug*(1-cosh(16._rprec*pi*zuvp(jz))/               &
        !     !                      2._rprec/(sinh(8._rprec*pi*lz_tot))**2)
        !     ! end if
        ! end do
    end if
else
    do jz = 1, nz
        ustar_ratio = ug*vonk/log(z_inv/zo1)
        arg2 = zuvp(jz)/(zo1/z_i)
        arg = ustar_ratio*(1._rprec/vonk)*log(arg2) !-1./(2.*vonk*z_i*z_i)*z*z
        if (coriolis_forcing) then
            ubar(jz) = ug
            vbar(jz) = vg
            wbar(jz) = 0._rprec
        else
            ubar(jz) = arg
            vbar(jz) = 0._rprec
            wbar(jz) = 0._rprec
        end if

        if (zuvp(jz) .gt. z_inv/z_i) then
            ubar(jz) = ug
        end if
    end do
end if 

! TOMAS CHOR added for stress-free IC
if ((coriol==0.0_rprec) .or. (u_star==0.0_rprec)) then
    ubar(:) = 0.0_rprec
    vbar(:) = 0.0_rprec
    wbar(:) = 0.0_rprec
end if

rms = 1._rprec
bv_freq = (g/T_scale*inv_strength)**(1._rprec/2._rprec)
do jz = 1, nz-1
    jz_abs = coordz*(nz - 1) + jz
    z = (coordz*(nz-1) + jz - 0.5_rprec)*dz*z_i
    seed = -80 - jz_abs !--trying to make consistent init for MPI
    do jy = 1, nynpy
    do jx = 1, nxt
        !...Ran3 returns uniform RV between 0 and 1. (sigma_rv=0.289)
        !...Taking std dev of vel as 1 at all heights
        if (z .le. z_inv .or. z .ge. (lz_tot*z_i-ug_dim/bv_freq)) then
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            u(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + ubar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            v(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + vbar(jz)
            noise = rms/.289_rprec*(ran3(seed) - .5_rprec)
            w(jx, jy, jz) = noise*(1._rprec - z/z_i)*w_star/u_scale + wbar(jz)
            !!!noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            !!!theta(jx, jy, jz) = (theta_mean + 10._rprec*noise*(1 - z/z_i)*T_star)/T_scale
            !!!noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            !!!q(jx, jy, jz) = q_mix + 50._rprec*noise*(1 - z/z_i)*q_star
        else
            u(jx, jy, jz) = ubar(jz)
            v(jx, jy, jz) = vbar(jz)
            w(jx, jy, jz) = wbar(jz)
            ! if ((z .gt. z_turb_limit) .and. (z .le. z_inv)) then
            !     theta(jx, jy, jz) = theta_mean/T_scale
            ! else
            !!!    theta(jx, jy, jz) = (theta_mean + (z - z_inv)*inv_strength)/T_scale
            !!!    ! theta(jx,jy,jz) = (theta_mean + (z-z_turb_limit)*inv_strength)/T_scale
            !!!! end if
            !!!q(jx, jy, jz) = q_mix
        end if
        
        if (z .le. z_inv) then
            noise = rms/0.289_rprec*(ran3(seed) - .5_rprec)
            theta(jx, jy, jz) = (theta_mean + 10._rprec*noise*(1 - z/z_i)*T_star)/T_scale
        else
            theta(jx, jy, jz) = (theta_mean + (z - z_inv)*inv_strength)/T_scale
        end if
    end do
    end do
end do
! ---
if ( coordz == 0 ) then
    w(1:nxt, 1:nynpy, 1) = 0._rprec
end if
if ( coordz == npz - 1 ) then
    w(1:nxt, 1:nynpy, nz) = 0._rprec
endif

if ( coordz == npz - 1 ) then
    if ( ubc_mom == 'wall' ) then
        u(1:nxt, 1:nynpy, nz) = BOGUS
        v(1:nxt, 1:nynpy, nz) = BOGUS
        theta(1:nxt, 1:nynpy, nz) = BOGUS
    else if ( ubc_mom == 'stress free' ) then
        u(1:nxt, 1:nynpy, nz) = u(1:nxt, 1:nynpy, nz - 1)
        v(1:nxt, 1:nynpy, nz) = v(1:nxt, 1:nynpy, nz - 1)
        theta(1:nxt, 1:nynpy, nz) = theta(1:nxt, 1:nynpy, nz - 1) + dTdz_top/T_scale*z_i*dz
    end if
end if
! ---
end subroutine ic_scal
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
function ran3(idum)
!-----------------------------------------------------------------------
!   random number generator
!-----------------------------------------------------------------------
use types, only:rprec
implicit none
! ---
integer(kind=4) :: idum
real(kind=rprec)    :: ran3
! ---
integer, parameter :: mbig  = 1000000000
integer, parameter :: mseed = 161803398
integer, parameter :: mz    = 0
real(kind=rprec), parameter :: fac = 1.d0/mbig
!parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
! ---
integer :: i, iff, ii, inext, inextp, k
integer :: mj, mk, ma(55)
save iff, inext, inextp, ma
data iff /0/
! ---  
if( idum .lt. 0 .or. iff .eq. 0 ) then
    iff = 1
    mj = mseed - iabs(idum)
    mj = mod(mj, mbig)
    ma(55) = mj
    mk = 1
    do i = 1, 54
        ii = mod(21*i, 55)
        ma(ii) = mk
        mk = mj - mk
        if(mk .lt. mz) mk = mk + mbig
        mj = ma(ii)
    end do
    
    do k = 1, 4
        do i = 1, 55
            ma(i) = ma(i) - ma(1+mod(i+30,55))
            if(ma(i) .lt. mz) ma(i) = ma(i) + mbig
        end do
    end do
    inext = 0
    inextp = 31
    idum = 1
end if
  
inext = inext + 1
if(inext .eq. 56) inext = 1
inextp = inextp + 1
if(inextp .eq. 56) inextp = 1
mj = ma(inext) - ma(inextp)
if(mj .lt. mz) mj = mj + mbig
ma(inext) = mj
ran3 = mj * fac
! ---
return
end function ran3
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use sgsmodule
use io, only:fcumulative_time
implicit none
! ---
character(80) :: fname
logical :: exst
! ---
real(rprec), dimension(ldx, nyt, nzt) :: u_tot
real(rprec), dimension(ldx, nyt, nzt) :: v_tot
real(rprec), dimension(ldx, nyt, nzt) :: w_tot
real(rprec), dimension(ldx, nyt, nzt) :: RHSx_tot
real(rprec), dimension(ldx, nyt, nzt) :: RHSy_tot
real(rprec), dimension(ldx, nyt, nzt) :: RHSz_tot
real(rprec), dimension(ldx, nyt, nzt) :: Cs_opt2_tot
real(rprec), dimension(ldx, nyt, nzt) :: F_LM_tot
real(rprec), dimension(ldx, nyt, nzt) :: F_MM_tot
real(rprec), dimension(ldx, nyt, nzt) :: F_QN_tot
real(rprec), dimension(ldx, nyt, nzt) :: F_NN_tot
! ---
if ( rank == 0 ) then
    inquire (file=fcumulative_time, exist=exst)
    if (exst) then
        open(1, file=fcumulative_time)
        read(1, *) nums
        close(1)
    else
        write(*, *) "Restart file NOT found"
        nums = 0
    end if

    write (fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    inquire (file=fname, exist=exst)
    if (exst) then
        open(1, file=fname, form='unformatted')
        print *, 'Reading initial velocity from file'
        read(1) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt),      &
                w_tot(:, :, 1:nzt), RHSx_tot(:, :, 1:nzt),   &
                RHSy_tot(:, :, 1:nzt), RHSz_tot(:, :, 1:nzt), Cs_opt2_tot, &
                F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
        close (1)
    else 
        print *, 'Cannot read initial velocity from file. Stop.'
        stop
    end if
end if
! ---
call scatter_data(u_tot, u)
call scatter_data(v_tot, v)
call scatter_data(w_tot, w)
call scatter_data(RHSx_tot, RHSx)
call scatter_data(RHSy_tot, RHSy)
call scatter_data(RHSz_tot, RHSz)
call scatter_data2(Cs_opt2_tot, Cs_opt2)
call scatter_data2(F_LM_tot, F_LM)
call scatter_data2(F_MM_tot, F_MM)
call scatter_data2(F_QN_tot, F_QN)
call scatter_data2(F_NN_tot, F_NN)
! ---
end subroutine read_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_vel_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use sgsmodule
use io, only:fcumulative_time
implicit none
! ---
character(80) :: fname
logical :: exst
! ---
real, dimension(ldx, nyt, nzt) :: u_tot
real, dimension(ldx, nyt, nzt) :: v_tot
real, dimension(ldx, nyt, nzt) :: w_tot

! ---
if ( rank == 0 ) then
    inquire (file=fcumulative_time, exist=exst)
    if (exst) then
        open(1, file=fcumulative_time)
        read(1, *) nums
        close(1)
    else
        write(*, *) "Restart file NOT found"
        nums = 0
    end if

    ! write (fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    write (fname, '(a)') '../readin/vel_tt.out'
    inquire (file=fname, exist=exst)
    if (exst) then
        open(1, file=fname, form='unformatted')
        print *, 'Reading initial velocity from file'
        read(1) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt),      &
                w_tot(:, :, 1:nzt)
        close (1)
    else 
        print *, 'Cannot read initial velocity from file. Stop.'
        stop
    end if
end if
! ---
call scatter_data3(u_tot, u)
call scatter_data3(v_tot, v)
call scatter_data3(w_tot, w)
! ---
end subroutine read_vel_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
!subroutine read_flow_field()
!!-----------------------------------------------------------------------
!!   random number generator
!!-----------------------------------------------------------------------
!use types, only:rprec
!use param
!use sim_param
!use sgsmodule
!use io, only:fcumulative_time
!implicit none
!! ---
!real(rprec), dimension(:, :, :), allocatable :: u_tot, v_tot, w_tot
!real(rprec), dimension(:, :, :), allocatable :: RHSx_tot, RHSy_tot, RHSz_tot
!real(rprec), dimension(:, :, :), allocatable :: Cs_opt2_tot
!real(rprec), dimension(:, :, :), allocatable :: F_LM_tot
!real(rprec), dimension(:, :, :), allocatable :: F_MM_tot
!real(rprec), dimension(:, :, :), allocatable :: F_QN_tot
!real(rprec), dimension(:, :, :), allocatable :: F_NN_tot
!
!integer, dimension(npz) :: sendcounts, recvcounts, displs
!character(80) :: fname
!logical :: exst
!! ---
!if ( rank == 0 ) then
!    inquire (file=fcumulative_time, exist=exst)
!    if (exst) then
!        open(1, file=fcumulative_time)
!        read(1, *) nums
!        close(1)
!    else
!        write(*, *) "Restart file NOT found"
!        nums = 0
!    end if
!
!    !-- allocate variable
!    allocate (u_tot(ldx, nyt, nzt),       v_tot(ldx, nyt, nzt),     &
!              w_tot(ldx, nyt, nzt),       RHSx_tot(ldx, nyt, nzt),  &
!              RHSy_tot(ldx, nyt, nzt),    RHSz_tot(ldx, nyt, nzt),  &
!              Cs_opt2_tot(ldx, nyt, nzt), F_LM_tot(ldx, nyt, nzt),  &
!              F_MM_tot(ldx, nyt, nzt),    F_QN_tot(ldx, nyt, nzt),  &
!              F_NN_tot(ldx, nyt, nzt))
!    !-- read
!    write(fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
!    inquire (file=fname, exist=exst)
!    if ( exst ) then
!        open(1, file=fname, form='unformatted')
!        write(*,*) 'Reading initial velocity from file'
!        read(1) u_tot, v_tot, w_tot, RHSx_tot, RHSy_tot, RHSz_tot,      &
!                Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
!        close(1)
!    else
!        write(*,*) 'Cannot read initial velocity from file. Stop.'
!        stop
!    end if
!end if
!!- scatter the data
!sendcounts = size(u(:, :, 1:nz-1))
!recvcounts = size(u(:, :, 1:nz))
!displs = coord_of_rank * sendcounts
!call mpi_scatterv(u_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      u(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(v_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      v(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(w_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      w(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(RHSx_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      RHSx(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(RHSy_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      RHSy(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(RHSz_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      RHSz(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(Cs_opt2_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      Cs_opt2(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(F_LM_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      F_LM(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(F_MM_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      F_MM(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(F_QN_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      F_QN(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!call mpi_scatterv(F_NN_tot(1, 1, 1), recvcounts, displs, MPI_RPREC, &
!                      F_NN(1, 1, 1), recvcounts, MPI_RPREC, rank_of_coord(0), comm, ierr)
!
!! ---
!if ( rank == 0 ) then
!    deallocate(u_tot, v_tot, w_tot, RHSx_tot, RHSy_tot, RHSz_tot)
!    deallocate(Cs_opt2_tot, F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot)
!end if
!
!end subroutine read_flow_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine ekman_layer(ubar, vbar)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
use param
use stokes_drift, only: U_stokes, wavenm_w
implicit none
! ---
integer :: jz
real(kind=rprec):: z, nu
real(kind=rprec), dimension(nz), intent(inout) :: ubar, vbar
complex(kind=rprec) :: gamma, Ubarcomp
! ---
nu = nu_ek/(u_scale*z_i)

do jz = 1, nz
    z = zuvp(jz)
    
    gamma = cmplx(0._rprec, 1._rprec)*coriol*U_stokes/            &
            (4._rprec*wavenm_w**2*nu - cmplx(0._rprec, 1._rprec)*coriol)

    Ubarcomp = cmplx(1._rprec, -1._rprec) / sqrt(2._rprec*coriol*nu) &
             * (u_star - 2._rprec*wavenm_w*nu*gamma &
             * cmplx(cos(agl_Ust), sin(agl_Ust)))*exp(cmplx(1._rprec, 1._rprec) &
             / sqrt(2._rprec)*sqrt(coriol/nu)*(-z)) &
             + cmplx(cos(agl_Ust), sin(agl_Ust))*gamma*exp(2._rprec*wavenm_w*(-z))
    !BC END Changed by Bicheng Chen for direction of Stokes drift
    ubar(jz) = real(Ubarcomp) + ug
    vbar(jz) = aimag(Ubarcomp) + vg
end do
! ---
end subroutine ekman_layer
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_profile(ubar, vbar, z_inv)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
use param
use stokes_drift, only: U_stokes, wavenm_w
implicit none
! ---
integer :: jz
real(kind=rprec):: z, nu, Kvis
real(kind=rprec), intent(in) :: z_inv
real(kind=rprec), dimension(nz), intent(inout) :: ubar, vbar
complex(kind=rprec) :: gamma, Ubarcomp
! ---
nu = nu_ek/(u_scale*z_i)
Kvis = 1.e-4_rprec

do jz = 1, nz
    z = zuvp(jz)
    
    if ( z .le. z_inv/z_i ) then
        gamma = cmplx(0._rprec, 1._rprec)*coriol*U_stokes/            &
                (4._rprec*wavenm_w**2*nu - cmplx(0._rprec, 1._rprec)*coriol)

        Ubarcomp = cmplx(1._rprec, -1._rprec) / sqrt(2._rprec*coriol*nu) &
                    * (u_star - 2._rprec*wavenm_w*nu*gamma &
                    * cmplx(cos(agl_Ust), sin(agl_Ust)))*exp(cmplx(1._rprec, 1._rprec) &
                    / sqrt(2._rprec)*sqrt(coriol/nu)*(-z)) &
                    + cmplx(cos(agl_Ust), sin(agl_Ust))*gamma*exp(2._rprec*wavenm_w*(-z))
        !BC END Changed by Bicheng Chen for direction of Stokes drift
        ubar(jz) = real(Ubarcomp) + ug
        vbar(jz) = aimag(Ubarcomp) + vg
    else
        gamma = sqrt(freq_coriolis/2._rprec/Kvis)
        ubar(jz) = ug*(1-exp(-gamma*z*z_i)*cos(gamma*z*z_i))
        vbar(jz) = ug*exp(-gamma*z*z_i)*sin(gamma*z*z_i)
    end if
end do
! ---
end subroutine init_profile
