!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init()
!--------------------------------------------------------------------!
!   initialize everything to make main program clean                                           
!--------------------------------------------------------------------!    
use param, only: theta_flag, pcon_flag, model, ubc, sgs, comm, ierr
use sim_param, only: u, v, w
use io, only: screen_display
use intermediate
use bottombc, only: patches
use fft
use test_filtermodule
use topbc, only: sponge, setsponge
implicit none
! ---
call init_namelist()            ! read all pararmeters from namelist 
! ---
call init_parameter()           ! domain variables initialization
! ---
call init_nondimensional()      ! normalize the variable
! ---
call allocate_output_variable() ! output variables initialization
call allocate_flow_variable()   ! flow variables initialization
! ---
call init_vel_field()
call init_temperature_field()
call init_concentration_field()
call data_exchange()            ! exchange data between neighboring zones
! ---
call patches()                  ! Initialize surface physical conditions
! --- formulate the fft plans
call init_fft_plan()            ! create fft plans and initialize the kx, ky arrays
! --- 
!call openfiles()                ! open output files
! --- initialize test filter
if ( sgs ) then
    call init_test_filter(2._rprec * filter_size, G_test)

    if (model == 3 .or. model == 5) then  !--scale dependent dynamic
        call init_test_filter(4._rprec * filter_size, G_test_test)
    end if
end if
! --- define the upper boundary condition (sponge damping layer)
if (ubc == 1) then
    call setsponge()
else
    sponge = 0._rprec
end if
! ---
call init_lad()
! --- 
call filter_data(u)
call filter_data(v)
call filter_data(w)
! ---
! call screen_display()           ! echo the simulation parameters
! ---
end subroutine init
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init_namelist()
!--------------------------------------------------------------------!
!     Read all the pararmeters in namelist                                        
!--------------------------------------------------------------------!
use param
implicit none
! ---
open(fid_param, file=fn_param, form='formatted', status='old')
read(fid_param, nml=output_control)
if (flagVSpl) then
    allocate (xVSplS(nVSpl), xVSplE(nVSpl), yVSplS(nVSpl), yVSplE(nVSpl), &
              zVSplS(nVSpl), zVSplE(nVSpl))
end if
if (flagCSpl) then
    allocate (xCSplS(nVSpl), xCSplE(nVSpl), yCSplS(nVSpl), yCSplE(nVSpl), &
              zCSplS(nVSpl), zCSplE(nVSpl))
end if
  
read (fid_param, nml=output_sample)
read (fid_param, nml=domain_param)
read (fid_param, nml=time_param)
read (fid_param, nml=flow_param)
read (fid_param, nml=canopy_param)
read (fid_param, nml=temperature_param)

read (fid_param, nml=ocean_param)  
read (fid_param, nml=con_param)

if (pcon_flag) then
    info_con%n_con => n_con
    if (settling) then
        allocate (vel_settling(n_con))
        allocate (ratio_dens(n_con))
        allocate (vol_spec(n_con))
        read(fid_param, nml=particle_param)
        info_con%vel_settling => vel_settling
        info_con%ratio_dens => ratio_dens
        info_con%vol_spec => vol_spec
    end if 
end if
close (fid_param)
! ---
end subroutine init_namelist
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_parameter()
!-----------------------------------------------------------------------
!   domain variables initialization (Bicheng Chen 07/02/2015)  
!-----------------------------------------------------------------------
use param
use scalars_module, only:hetero_count_out
implicit none
! ---
integer :: ind
! --- 
nz   = (nzt - 1)/npz + 1

nxt2  = 3*nxt/2
nyt2  = 3*nyt/2
lhx  = nxt/2 + 1
lhy  = nyt/2 + 1
ldx  = 2*lhx
ldy  = 2*lhy
lhx2 = nxt2/2 + 1
ldx2 = 2*lhx2
lhy2 = nyt2/2 + 1
ldy2 = 2*lhy2

nxnpy  = nxt / npy
nynpy  = nyt / npy
nx2npy = nxt2/ npy
ny2npy = nyt2/ npy
nxhnpy = nxnpy/2 + 1
nxh2npy= nx2npy/2 + 1

dx   = lx_tot / nxt
dy   = ly_tot / nyt
dz   = lz_tot / (nzt - 1)

ly   = ly_tot / npy 
lz   = lz_tot / npz

! --- initialize other time variables
hetero_count_out = p_count
! ---
! if (flag_restart) then
!     info_time%flag_restart = flag_restart
!     info_time%nt_restart = nt_restart
! end if
! ---
end subroutine init_parameter
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_nondimensional()
!-----------------------------------------------------------------------
!   normalize the variables
!-----------------------------------------------------------------------
use types, only:rprec
use param
use stokes_drift
use intermediate
implicit none
! --- 
ly = ly / z_i
lz = lz / z_i

lx_tot = lx_tot/z_i
ly_tot = ly_tot/z_i
lz_tot = lz_tot/z_i

dx = dx/z_i
dy = dy/z_i
dz = dz/z_i

call init_meshgrid()
! --- 
dt = dt * u_scale / z_i
u_star = u_star / u_scale
! --- canopy normalization
if ( use_res_canopy ) then
    hc  = hc  / z_i
    hgap = hgap / z_i
    lad = lad * z_i
end if
! --- Coriolis effect
if (coriolis_forcing) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
end if
! --- Mean pressure gradient force
if (use_mean_p_force) then
    if (mean_p_force == float_missing) then
        mean_p_force = 1._rprec / lz_tot
    else
        mean_p_force = mean_p_force / u_scale**2 * z_i
    end if
else 
    mean_p_force = 0._rprec
end if

! ---
if (ocean_flag .and. theta_flag) then
    alpha_w = alpha_w * T_scale     ! normalize thermal expansion rate
end if

! --- Ocean flag
call init_stokes_drift

if (ocean_flag) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
    !- if use dynamical stress, read friction velocity and its angle
    if (flag_dynStress) then
        inquire (iolength=len_ustar) ustar_dyn, agl_stress
        open (unit=fid_ustar, file=fn_ustar, access='direct', recl=len_ustar)
        read (unit=fid_ustar, rec=1) ustar_dyn, agl_stress
        close (fid_ustar)
        ustar_dyn = ustar_dyn / u_scale
    end if 
    rad_stress = agl_stress/180._rprec*pi
    
    ! >>Cross flow over all domain
    if (flag_crossVel) then
        u_cross = udim_cross / u_scale
        v_cross = vdim_cross / u_scale
    end if
end if
! ---
end subroutine init_nondimensional
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_meshgrid()
!-----------------------------------------------------------------------
!    
!-----------------------------------------------------------------------
use param
implicit none
! ---
integer :: jx, jy, jz
! ---
allocate(x(ldx), y(nynpy), zuvp(0:nz), zw(0:nz))
! --- 
do jz = 0, nz
    zw(jz)   = (coordz*(nz-1) + jz - 1) * dz
    zuvp(jz) = (coordz*(nz-1) + jz - 0.5_rprec) * dz
end do

do jx = 1, ldx 
    x(jx) = (jx - 1) * dx
end do

do jy = 1, nynpy
    y(jy) = (coordy*nynpy + jy - 1) * dy
end do
! ---
end subroutine init_meshgrid
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_output_variable()
!-----------------------------------------------------------------------
!   output variables initialization (Bicheng Chen 06/13/2016)
!-----------------------------------------------------------------------
use io
use param !, only:nxt, nyt, nz, nzt, use_avgslice, average_dim_num, USE_MPI
use intermediate !, only:average_dim_select_flag, dim1_size, dim2_size, &
                 !     dim1_global, dim2_global, avg_out
implicit none
! ---
average_dim_select_flag = 1 - average_dim_num/2
dim1_size   = average_dim_select_flag*(nxt - nz + 1) + nz - 1
dim2_size   = average_dim_select_flag*(nz - 2) + 1
dim1_global = average_dim_select_flag*(nxt - nzt + 1) + nzt - 1
dim2_global = average_dim_select_flag*(nzt - 2) + 1
  
if (use_avgslice) then
    if (average_dim_num == 1) then
        allocate (avg_out(1:nxt, 1:nzt - 1))
    else if (average_dim_num == 2) then
        allocate (ax_1d(nz - 1))
        allocate (avg_out_1d(nzt - 1))
    end if
    allocate (ap(nxt, nz - 1), au(nxt, nz - 1), av(nxt, nz - 1), aw(nxt, nz - 1),       &
              p2(nxt, nz - 1), u2(nxt, nz - 1), v2(nxt, nz - 1), w2(nxt, nz - 1),       &
              auw(nxt, nz - 1), avw(nxt, nz - 1), acs(nxt, nz - 1), adudz(nxt, nz - 1), &
              advdz(nxt, nz - 1), aCs_Ssim(nxt, nz - 1), abeta_sgs(nxt, nz - 1),       &
              abetaclip_sgs(nxt, nz - 1), atxx(nxt, nz - 1), atxz(nxt, nz - 1),        &
              atyy(nxt, nz - 1), atyz(nxt, nz - 1), atzz(nxt, nz - 1),                 &
              u3(nxt, nz - 1), v3(nxt, nz - 1), w3(nxt, nz - 1),                       &
              adudt(nxt, nz - 1), advdt(nxt, nz - 1), adwdt(nxt, nz - 1))
    allocate (awu2(nxt, nz - 1), awv2(nxt, nz - 1), awT2(nxt, nz - 1),     &
              auv(nxt, nz - 1), awp(nxt, nz - 1), autau13(nxt, nz - 1),    &
              avtau23(nxt, nz - 1), adissip(nxt, nz - 1))
    
    if (theta_flag) then 
        allocate (atheta(nxt, nz - 1), t2(nxt, nz - 1), q2(nxt, nz - 1),       &
                  asgs_t3(nxt, nz - 1), awt(nxt, nz - 1), adTdz(nxt, nz - 1),  &
                  anu_t(nxt, nz - 1), t3(nxt, nz - 1), var_t(nxt, nz - 1))
    end if
    
    if (PCon_FLAG) then 
        allocate (aPCon(nxt, nz - 1), PCon2(nxt, nz - 1), asgs_PCon3(nxt, nz - 1),     &
                  awPCon(nxt, nz - 1), adPCondz(nxt, nz - 1), var_PCon(nxt, nz - 1),   &
                  aKc_t(nxt, nz - 1), aCs2Sc(nxt, nz - 1), asgs_PCon1(nxt, nz - 1),    & 
                  auPCon(nxt, nz - 1))
    end if
end if
  
if ( rank == 0 ) then
    allocate (ubar_avg(1, nzt - 1), vbar_avg(1, nzt - 1),             &
              thetabar_avg(1, nzt - 1), Cs2bar_avg(1, nzt - 1),       &
              Nutbar_avg(1, nzt - 1), ubar_tot(1, nzt - 1),           &
              vbar_tot(1, nzt - 1), thetabar_tot(1, nzt - 1),         &
              Cs2bar_tot(1, nzt - 1), Nutbar_tot(1, nzt - 1))
    allocate (upr2bar_avg(3, nzt - 1), stressbar_avg(3, nzt - 1),     &
              upr2bar_tot(3, nzt - 1), stressbar_tot(3, nzt - 1))
    allocate (Eozbar_avg(1, nxt/2, nzt - 1), Eozbar_tot(1, nxt/2, nzt - 1))
end if

! --- io module
!jx_pls = 1; jx_ple = nxt; width = nyt/2 - 1
!jy_pls = nyt/2 - width; jy_ple = nyt/2 + width + 1
!allocate(mean_u(jx_pls:jx_ple, jy_pls:jy_ple, nz),  &
!         mean_v(jx_pls:jx_ple, jy_pls:jy_ple, nz),  &
!         mean_w(jx_pls:jx_ple, jy_pls:jy_ple, nz),  &
!         mean_u2(jx_pls:jx_ple, jy_pls:jy_ple, nz), &
!         mean_v2(jx_pls:jx_ple, jy_pls:jy_ple, nz), &
!         mean_w2(jx_pls:jx_ple, jy_pls:jy_ple, nz))
!
!! ---
!mean_u = 0._rprec; mean_u2 = 0._rprec
!mean_v = 0._rprec; mean_v2 = 0._rprec
!mean_w = 0._rprec; mean_w2 = 0._rprec

! ---
end subroutine allocate_output_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_flow_variable()
!-----------------------------------------------------------------------
!   flow variables initialization (Bicheng Chen 06/07/2016)
!! >>Add dynamical wind stress and dynamical Stokes drift
!! >>(Bicheng Chen 10/20/2016)
!-----------------------------------------------------------------------
use param
use sim_param
use fft
use test_filtermodule, only:G_test, G_test_test
use sgsmodule
use bottombc, only:ustar_avg, zo, z_os, patch, d0,          &
                   T_s, q_s, phi_m, psi_m, phi_h, psi_h,    &
                   zo_PCon, PCon_s
use topbc, only:sponge
use intermediate
use scalars_module
implicit none
! ---
integer :: jx, jy, ind

integer, dimension(npz) :: sendcounts, recvcounts, displs

real(rprec), dimension(:, :, :, :), allocatable :: k_ttlCR_tot
character(20), dimension(:), allocatable :: unt_bg      ! initial unit of background concentration
character(20) :: unt_to  
real(rprec), dimension(:, :), allocatable :: con_bg
real(rprec), dimension(:, :, :, :), allocatable :: con_bg_3d
real(rprec), dimension(:, :), allocatable :: k_bg       ! chemical reaction rate of specified species
integer, parameter :: fid_bg = 822
namelist/conprof/  unt_bg, unt_to, con_bg
namelist/chemrate/ k_bg
! ---
allocate (u(ldx, nynpy, 0:nz), v(ldx, nynpy, 0:nz), w(ldx, nynpy, 0:nz))
allocate (ke(ldx, nynpy, 0:nz), ke_temp(ldx, nynpy, 0:nz), &
          dragx(ldx, nynpy, 0:nz), dragy(ldx, nynpy, 0:nz), dragz(ldx, nynpy, 0:nz))
allocate (dudx(ldx, nynpy, 0:nz), dudy(ldx, nynpy, 0:nz), dudz(ldx, nynpy, 0:nz),        &
          dvdx(ldx, nynpy, 0:nz), dvdy(ldx, nynpy, 0:nz), dvdz(ldx, nynpy, 0:nz),        &
          dwdx(ldx, nynpy, 0:nz), dwdy(ldx, nynpy, 0:nz), dwdz(ldx, nynpy, 0:nz),        &
          RHSx(ldx, nynpy, 0:nz), RHSy(ldx, nynpy, 0:nz), RHSz(ldx, nynpy, 0:nz),        &
          RHSx_f(ldx, nynpy, 0:nz), RHSy_f(ldx, nynpy, 0:nz), RHSz_f(ldx, nynpy, 0:nz))
allocate (dudt(ldx, nynpy, 0:nz), dvdt(ldx, nynpy, 0:nz), dwdt(ldx, nynpy, 0:nz))
allocate (dkedx(ldx, nynpy, 0:nz), dkedy(ldx, nynpy, 0:nz), dkedz(ldx, nynpy, 0:nz))
allocate (txx(ldx, nynpy, 0:nz), txy(ldx, nynpy, 0:nz), txz(ldx, nynpy, 0:nz),           &
          tyy(ldx, nynpy, 0:nz), tyz(ldx, nynpy, 0:nz), tzz(ldx, nynpy, 0:nz))
allocate (divtx(ldx, nynpy, 0:nz), divty(ldx, nynpy, 0:nz), divtz(ldx, nynpy, 0:nz))
allocate (p(ldx, nynpy, 0:nz))
allocate (dpdx(ldx, nynpy, nz), dpdy(ldx, nynpy, nz), dpdz(ldx, nynpy, nz))

allocate (kx(lhx), ky(lhx))
allocate (kx_2d(lhx, nyt), ky_2d(lhx, nyt), k2(lhx, nyt))
allocate (kx_2d_mpi(nxhnpy, nyt), ky_2d_mpi(nxhnpy, nyt), k2_mpi(nxhnpy, nyt))
allocate (dealias(nxhnpy, nyt))
allocate (G_test(nxhnpy, nyt), G_test_test(nxhnpy, nyt))

allocate (Nu_t(ldx, nynpy, 0:nz), dissip(ldx, nynpy, 0:nz), magS(ldx, nynpy, 0:nz),     &
          u_lag(ldx, nynpy, 0:nz),v_lag(ldx, nynpy, 0:nz), w_lag(ldx, nynpy, 0:nz))
allocate (F_LM(ldx, nynpy, nz), F_MM(ldx, nynpy, nz), F_QN(ldx, nynpy, nz),             &
          F_NN(ldx, nynpy, nz), Beta(ldx, nynpy, nz), Betaclip(ldx, nynpy, nz))
allocate (Beta_avg(nz), Betaclip_avg(nz))
allocate (Cs_opt2(ldx, nynpy, nz), Cs_opt2_avg(ldx, nynpy, nz), Cs_Ssim(ldx, nynpy, nz), &
          epsilon_lag(ldx, nynpy, nz), epsilon_lag2(ldx, nynpy, nz),                  &
          xlag(ldx, nynpy, nz), ylag(ldx, nynpy, nz), zlag(ldx, nynpy, nz))
allocate (F_KX(nxt, nynpy, nz), F_XX(nxt, nynpy, nz), F_KX2(nxt, nynpy, nz), F_XX2(nxt, nynpy, nz))

allocate (ustar_avg(nxt, nynpy), zo(nxt, nynpy), z_os(nxt, nynpy), patch(nxt, nynpy), d0(nxt, nynpy))
allocate (sponge(0:nz))

allocate (cross_x(ldx, nynpy, nz), cross_y(ldx, nynpy, nz), cross_z(ldx, nynpy, nz))
allocate (vort1(ldx, nynpy, 0:nz), vort2(ldx, nynpy, 0:nz), vort3(ldx, nynpy, 0:nz))
allocate (beta_scal(ldx, nynpy, 0:nz))

allocate (cross_x_big(nxt2, ny2npy, nz), cross_y_big(nxt2, ny2npy, nz), cross_z_big(nxt2, ny2npy, nz))
allocate (u_big(nxt2, ny2npy, 0:nz), v_big(nxt2, ny2npy, 0:nz), w_big(nxt2, ny2npy, 0:nz))
allocate (vort1_big(nxt2, ny2npy, 0:nz), vort2_big(nxt2, ny2npy, 0:nz), vort3_big(nxt2, ny2npy, 0:nz))

! --- press_stag_array
allocate (rH_x(ldx, nynpy, 0:nz), rH_y(ldx, nynpy, 0:nz), rH_z(ldx, nynpy, 0:nz))
allocate (H_x(nxhnpy, nyt, 0:nz), H_y(nxhnpy, nyt, 0:nz), H_z(nxhnpy, nyt, 0:nz))
allocate (rtopw(ldx, nynpy), rbottomw(ldx, nynpy))
allocate (topw(nxhnpy, nyt), bottomw(nxhnpy, nyt))

! --- SGS model
if ( sgs ) then
    if (model == 3) then
        allocate (L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                        &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                        &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                        &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy), S_bar(ldx, nynpy),     &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),            &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),            &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),      &
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),      &
                  S_hat(ldx, nynpy),       &
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),            &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),            &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),      &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),      &
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),                  &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy), S(ldx, nynpy))
        allocate (beta_sd(nz))
        M11 => Q11; M12 => Q12; M13 => Q13; M22 => Q22; M23 => Q23; M33 => Q33
    else if (model == 4) then
        allocate (u_pr(nxt, nynpy), v_pr(nxt, nynpy), w_pr(nxt, nynpy))
        allocate (w_nod(ldx, nynpy), S(ldx, nynpy), tempos(ldx, nynpy),                 &    
                  L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                    &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                    &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                    &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy),                    &
                  M11(ldx, nynpy), M12(ldx, nynpy), M13(ldx, nynpy),                    &
                  M22(ldx, nynpy), M23(ldx, nynpy), M33(ldx, nynpy),                    &
                  N11(ldx, nynpy), N12(ldx, nynpy), N13(ldx, nynpy),                    &
                  N22(ldx, nynpy), N23(ldx, nynpy), N33(ldx, nynpy),                    & 
                  LM(ldx, nynpy), MM(ldx, nynpy), QN(ldx, nynpy), NN(ldx, nynpy),       &
                  Tn(ldx, nynpy), epsi(ldx, nynpy), dumfac(ldx, nynpy), S_bar(ldx, nynpy),  &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),        &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),        &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),  &
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),  &
                  S_hat(ldx, nynpy),    & 
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),        &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),        &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),  &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),  &
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),              &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy), fourbeta(ldx, nynpy))
        allocate (xp(ldx, nynpy, 0:nz), yp(ldx, nynpy, 0:nz), zp(ldx, nynpy, 0:nz),     &
                  u_temp(ldx, nynpy, 0:nz), v_temp(ldx, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2))
  
    else if (model == 5) then
        allocate (visc(ldx, nynpy, nz), Cs_opt2_2d(ldx, nynpy, nz), Cs_opt2_4d(ldx, nynpy, nz))
        allocate (S(ldx, nynpy), tempos(ldx, nynpy),                                        & 
                  L11(ldx, nynpy), L12(ldx, nynpy), L13(ldx, nynpy),                        &
                  L22(ldx, nynpy), L23(ldx, nynpy), L33(ldx, nynpy),                        &
                  Q11(ldx, nynpy), Q12(ldx, nynpy), Q13(ldx, nynpy),                        &
                  Q22(ldx, nynpy), Q23(ldx, nynpy), Q33(ldx, nynpy),                        &
                  M11(ldx, nynpy), M12(ldx, nynpy), M13(ldx, nynpy),                        &
                  M22(ldx, nynpy), M23(ldx, nynpy), M33(ldx, nynpy),                        &
                  N11(ldx, nynpy), N12(ldx, nynpy), N13(ldx, nynpy),                        &
                  N22(ldx, nynpy), N23(ldx, nynpy), N33(ldx, nynpy),                        & 
                  LM(ldx, nynpy), MM(ldx, nynpy), QN(ldx, nynpy), NN(ldx, nynpy),           &
                  Tn(ldx, nynpy), epsi(ldx, nynpy), dumfac(ldx, nynpy), S_bar(ldx, nynpy),  &
                  S11_bar(ldx, nynpy), S12_bar(ldx, nynpy), S13_bar(ldx, nynpy),            &
                  S22_bar(ldx, nynpy), S23_bar(ldx, nynpy), S33_bar(ldx, nynpy),            &
                  S_S11_bar(ldx, nynpy), S_S12_bar(ldx, nynpy), S_S13_bar(ldx, nynpy),      & 
                  S_S22_bar(ldx, nynpy), S_S23_bar(ldx, nynpy), S_S33_bar(ldx, nynpy),      &
                  S_hat(ldx, nynpy),        & 
                  S11_hat(ldx, nynpy), S12_hat(ldx, nynpy), S13_hat(ldx, nynpy),            &
                  S22_hat(ldx, nynpy), S23_hat(ldx, nynpy), S33_hat(ldx, nynpy),            &
                  S_S11_hat(ldx, nynpy), S_S12_hat(ldx, nynpy), S_S13_hat(ldx, nynpy),      &
                  S_S22_hat(ldx, nynpy), S_S23_hat(ldx, nynpy), S_S33_hat(ldx, nynpy),      &  
                  u_bar(ldx, nynpy), v_bar(ldx, nynpy), w_bar(ldx, nynpy),                  &
                  u_hat(ldx, nynpy), v_hat(ldx, nynpy), w_hat(ldx, nynpy))
        allocate (LMvert(nz), MMvert(nz), QNvert(nz), NNvert(nz))
        allocate (xp(ldx, nynpy, 0:nz), yp(ldx, nynpy, 0:nz), zp(ldx, nynpy, 0:nz),         &
                  u_temp(ldx, nynpy, 0:nz), v_temp(ldx, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2),     &
                  FF_QN(nxt + 2, nynpy + 2, nz + 2), FF_NN(nxt + 2, nynpy + 2, nz + 2),     &
                  Beta_t(nxt + 2, nynpy + 2, nz + 2)) 
    end if
end if

! --- Temperature field
allocate (phi_m(nxt, nynpy), psi_m(nxt, nynpy), phi_h(nxt, nynpy), psi_h(nxt, nynpy))
allocate (L(nxt, nynpy), wstar(nxt, nynpy))

if (theta_flag) then
    allocate (T_s(nxt, nynpy), q_s(nxt, nynpy))
    allocate (theta(ldx, nynpy, 0:nz), q(ldx, nynpy, 0:nz))
    allocate (Pr_(ldx, nynpy, 0:nz), dTdz(ldx, nynpy, 0:nz), dqdz(ldx, nynpy, 0:nz),         &
              RHS_Tf(ldx, nynpy, 0:nz), RHS_T(ldx, nynpy, 0:nz), RHS_qf(ldx, nynpy, 0:nz),   &
              RHS_q(ldx, nynpy, 0:nz), sgs_t3(ldx, nynpy, 0:nz), sgs_q3(ldx, nynpy, 0:nz))    
    allocate (T_s_filtered(ldx, nynpy))
    
    ! - SGS model
    if (model_sc == 4) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2))
    else if (model_sc == 5) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2), &
                  FF_KX2(nxt + 2, nynpy + 2, nz + 2), FF_XX2(nxt + 2, nynpy + 2, nz + 2))
    end if
    !- use dynamical kinematic heat flux (Bicheng Chen)
    if (flag_dynWT) then
        inquire (iolength=len_wt) wt_s
        open (unit=fid_wt, file=fn_wt, access='direct', recl=len_wt)
        read (unit=fid_wt, rec=1) wt_s
        close (fid_wt)
    end if
end if

! --- concentration variables initialization
allocate (beta_pcon(ldx, 0:nynpy+1, 0:nz))  ! allocate the variable always needed
if (pcon_flag) then
    npcon = info_con%n_con
    allocate (zo_PCon(nxt, nynpy), PCon_s(nxt, nynpy))
    allocate (P_surf_flux(nxt, nynpy), deposition(nxt, nynpy),  &
              P_surf_flux_dep(nxt, nynpy), Real_dep(nxt, nynpy))
    allocate (matrix_x(nxt, nxt), matrix_y(nynpy, nynpy))
    allocate (dvector_x(nxt), dvector_y(nynpy))

    if (settling) then
        allocate (settling_vel(npcon), densratio_pcon(npcon), v_pcon(npcon))
        settling_vel = info_con%vel_settling(1:npcon)/u_scale
        densratio_pcon = info_con%ratio_dens(1:npcon)
        v_pcon = info_con%vol_spec(1:npcon)*pcon_scale
    end if

    allocate (Kc_t(ldx, 0:nynpy+1, 0:nz), dPCondz(ldx, 0:nynpy+1, 0:nz),          &
              sgs_PCon3(ldx, 0:nynpy+1, 0:nz), res_PCon3(ldx, 0:nynpy+1, 0:nz),   &
              sgs_PCon1(ldx, 0:nynpy+1, 0:nz), Cs2Sc(ldx, 0:nynpy+1, 0:nz))
    
    !- allocate concentration
    allocate (PCon(ldx, 0:nynpy+1, 0:nz, npcon))
    allocate (RHS_PConf(ldx, 0:nynpy+1, 0:nz, npcon), RHS_PCon(ldx, 0:nynpy+1, 0:nz, npcon))

    call read_source_info()
end if

! ---
Cs_opt2_avg = 0._rprec
! ---
end subroutine allocate_flow_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_vel_field()
!-----------------------------------------------------------------------
!   initialize velocity field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
implicit none
! ---
if ( restart_vel ) then
    call read_flow_field()
else
    if (theta_flag) then
        call ic_scal()
    else
        call ic()
    end if
end if

!- add cross velocity to the all velocity field (Bicheng Chen 09/07/2016)
if (flag_crossVel .and. ocean_flag) then
    u(1:nxt, 1:nynpy, 0:nz) = u(1:nxt, 1:nynpy, 0:nz) + u_cross
    v(1:nxt, 1:nynpy, 0:nz) = v(1:nxt, 1:nynpy, 0:nz) + v_cross
end if

! ---
end subroutine init_vel_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_temperature_field()
!-----------------------------------------------------------------------
!   initialize temperature field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:RHS_T, sgs_t3
use bottombc, only:psi_m
implicit none
! ---
real(rprec), dimension(:, :, :), allocatable :: theta_tot
real(rprec), dimension(:, :, :), allocatable :: RHS_T_tot

character(80) :: fname
logical :: exst
! ---
if (theta_flag .and. restart_theta .and. restart_vel) then
    if ( rank == 0 ) then
        !-- allocate variable
        allocate (theta_tot(ldx, nyt, nzt), RHS_T_tot(ldx, nyt, nzt))
        !-- read
        write (fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        inquire (file=fname, exist=exst)
        if (exst) then
            open(1, file=fname, form='unformatted')
            print *, 'Reading initial temperature from file'
            read(1) theta_tot(:, :, 1:nzt), RHS_T_tot(:, :, 1:nzt), &
                    sgs_t3(:, :, 1), psi_m
            close(1)
        else
            print *, 'Cannot read initial temperature from file. Stop.'
            stop
        end if
    end if
    
    !- scatter the data
    call scatter_data(theta_tot, theta)
    call scatter_data(RHS_T_tot, RHS_T)

    if ( rank == 0 ) then
        deallocate(theta_tot, RHS_T_tot)
    end if

else if (theta_flag .and. restart_theta .and. (.not. restart_vel)) then !LV1
    print *, "Cannot initialize temperature field with initializing velocity."
    print *, "Stop"
    stop
else if ((.not. theta_flag) .and. restart_theta .and. restart_vel) then !LV1
    print *, "Cannot initialize temperature field with theta_flag."
    print *, "Stop"
    stop
else
    !- the initialization of temperature without file is in ic_scal
end if
! ---
end subroutine init_temperature_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_concentration_field()
!-----------------------------------------------------------------------
!   concentration variables initialization
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:deposition, Real_dep, RHS_PCon, Kc_t
use sgsmodule, only:F_KX, F_XX, F_KX2, F_XX2
implicit none
! --- initialize the variable always needed
real(rprec), dimension(:, :, :, :), allocatable :: PCon_tot
real(rprec), dimension(:, :, :, :), allocatable :: RHS_PCon_tot
real(rprec), dimension(:, :, :), allocatable :: Kc_t_tot
real(rprec), dimension(:, :, :), allocatable :: F_KX_tot
real(rprec), dimension(:, :, :), allocatable :: F_XX_tot
real(rprec), dimension(:, :, :), allocatable :: F_KX2_tot
real(rprec), dimension(:, :, :), allocatable :: F_XX2_tot
! ---
character(80) :: fname
logical :: exst
integer :: ipcon
! ---
if (pcon_flag .and. restart_pcon) then
    if ( rank == 0 ) then
        !-- allocate variable
        allocate (PCon_tot(ldx, nyt, nzt, npcon),     &
                  RHS_PCon_tot(ldx, nyt, nzt, npcon))
        allocate (Kc_t_tot(ldx, nyt, nzt))
        allocate (F_KX_tot(nxt, nyt, nzt),  F_XX_tot(nxt, nyt, nzt),  &
                  F_KX2_tot(nxt, nyt, nzt), F_XX2_tot(nxt, nyt, nzt))
        !-- read initial concentration from file
        write (fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
        inquire (file=fname, exist=exst)
        if (exst) then
            open(1, file=fname, form='unformatted')
            print *, 'Reading initial concentration from file'
            read(1) PCon_tot(:, :, 1:nzt, 1:npcon), &
                    RHS_PCon_tot(:, :, 1:nzt, 1:npcon)
            close(1)
        else
            print *, 'Cannot read initial concentration from file. Stop.'
            stop
        end if
        !-- read other initial variable releated to concentration from file
        write (fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
        inquire (file=fname, exist=exst)
        if (exst) then
            open(1, file=fname, form='unformatted')
            print *, 'Reading initial concentration from file'
            read(1) deposition(:, :), Real_dep(:, :), &
                    Kc_t_tot, F_KX_tot, F_XX_tot, F_KX2_tot, F_XX2_tot
            close(1)
        else
            print *, 'Cannot read initial other concentration variables from file. Stop.'
            stop
        end if
    end if
    !- scatter the data
    do ipcon = 1, npcon
        call scatter_data4(PCon_tot, PCon)
        call scatter_data4(RHS_PCon_tot, RHS_PCon)
    end do
    call scatter_data4(Kc_t_tot, Kc_t)
    call scatter_data2(F_KX_tot, F_KX)
    call scatter_data2(F_XX_tot, F_XX)
    call scatter_data2(F_KX2_tot, F_KX2)
    call scatter_data2(F_XX2_tot, F_XX2)

    if ( rank == 0 ) then
        deallocate (PCon_tot, RHS_PCon_tot)
        deallocate (Kc_t_tot, F_KX_tot, F_XX_tot, F_KX2_tot, F_XX2_tot)
    end if
else if ((.not. pcon_flag) .and. restart_pcon) then
    print *, "Cannot initialize concentration field without pcon_flag. Stop."
    stop
else if (pcon_flag .and. .not. restart_pcon) then
    !- Always initialize pollen concentration with zeros
    !- Chamecki - 08/01/2006
    !write(*,*) 'Pollen profile initialized with zeros'
    PCon(:, :, :, :) = 0._rprec
    if ( coordz == 0 ) then
        PCon(:, :, 0, :) = BOGUS
    end if
end if
! ---
end subroutine init_concentration_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_source_info
!-----------------------------------------------------------------------
!   initialize point sources and single sources (Bicheng Chen 05/22/2015)
!-----------------------------------------------------------------------
use param
use scalars_module
implicit none
! ---
character(80) :: fname
integer :: ipcon
! ---
if (pointsource) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        !print *, 'n_src ', n_src
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        con_src(ipcon)%n = n_src
        con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
        con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
        con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
        con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
    end do
  
    ! if (.not. single_point) then
    !     do ipcon = 1, npcon
    !         if (con_src(ipcon)%n /= 1) then
    !             print *, "the number of source is not 1 for non-single point sources."
    !             stop
    !         end if
    !     end do
    ! end if
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do

else if (flag_1Dsrc) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        read(fid_src, nml=source_alignment)
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        select case (src_align)
        case ('lateral')
            con_src(ipcon)%n = n_src
            con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
            ! con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
            ! ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
            ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
            ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        case ('streamwise')
            con_src(ipcon)%n = n_src
            ! con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
            con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            ! ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
            ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
            ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
            ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        case default
            write (*, *) 'invalid source alignment type'
            stop
        end select
    end do

    !-- calculate the specific cpu coord for point sources
    !-- and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do

else if (flag_2Dsrc) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write (fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        !print *, 'n_src ', n_src
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        con_src(ipcon)%n = n_src
        con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
    end do

    !-- calculate the specific cpu coord for point sources
    !-- and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
    end do
end if
! ---
end subroutine read_source_info
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------