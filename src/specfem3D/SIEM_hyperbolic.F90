!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

module siem_hyperbolic

! Hyperbolic ("artificially-slowed-gravity") solver for the perturbed
! gravitational potential, as an alternative to the elliptic SIEM Poisson solve.
!
! SIEM forms the Level-1 system  K phi = F  and inverts it every time step. Here
! the elliptic constraint is recast as a damped wave (telegraph) equation
!
!     M phi_tt + 2 kappa M phi_t + cg^2 K phi = cg^2 F
!
! and marched explicitly, sub-cycled inside each elastic step. Its steady state
! is exactly K phi = F, so the scheme relaxes towards the Poisson solution
! rather than solving for it; the error is O((v/cg)^2) = O(1/kg^2).
!
! Hirai et al. (2016, Phys. Rev. D 93, 083006) introduced the undamped form;
! Maeda et al. (2024, MNRAS 527, 471) added the damping that makes relaxation,
! rather than propagation, the intended behaviour.
!
! Structurally this mirrors the elastic and outer-core solvers rather than the
! SIEM CG solver: K is never stored or assembled, but applied element by element
! through the tensor-product GLL operator, exactly as compute_forces_* does for
! displacement. The mass matrix is diagonal and inverted by multiplication, and
! the damping term plays the role of the Stacey and Coriolis velocity terms.
!
! Two things are not matrix-free:
!  - the infinite-element layer, whose Zienkiewicz mapped-infinite shape
!    functions are not a tensor-product basis. It keeps its stored element
!    matrices, but it is a single layer of surface elements so the cost is
!    negligible.
!  - the load vector F, which still comes from compute_poisson_load3().

  use constants, only: CUSTOM_REAL,NDIM,myrank,IMAIN,ADD_TRINF, &
    NGLLX,NGLLY,NGLLZ,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  implicit none

  private

  public :: hyperbolic_prepare
  public :: hyperbolic_solve

  ! wavefield state on the Level-1 global degrees of freedom
  ! (the potential itself lives in pgrav1)
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: phi_dot,phi_ddot

  ! inverse of the assembled lumped mass matrix
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: minv1

  ! load vector at the previous elastic step, for source interpolation
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: load_prev

  ! work space
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: work_res,work_acc

  ! 1D GLL operators of the Level-1 basis
  ! hprime_gll(i,l) = derivative of basis l at point i, as in hprime_xx
  real(kind=CUSTOM_REAL),dimension(NGLLX_INF,NGLLX_INF) :: hprime_gll,hprimew_gll
  real(kind=CUSTOM_REAL),dimension(NGLLX_INF) :: wgll1d

  ! metric factors of the Laplacian, J * (contravariant metric), packed as
  ! (g11,g12,g13,g22,g23,g33) per Level-1 GLL point. This replaces the stored
  ! element stiffness matrices: 6 reals per point instead of NGLLCUBE_INF.
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: gfac_cm,gfac_oc,gfac_ic,gfac_trinf

  ! stride from Level-1 to Level-2 GLL points (2 for GLL3 inside GLL5)
  integer :: igll_stride

  ! artificial gravity speed squared and damping rate (non-dimensional)
  real(kind=CUSTOM_REAL) :: cg2,kappa

  ! sub-stepping
  real(kind=CUSTOM_REAL) :: dt_sub
  integer :: nsub

contains

!
!-------------------------------------------------------------------------------
!

  subroutine hyperbolic_prepare()

  use constants, only: HYPERBOLIC_KG,HYPERBOLIC_KAPPA, &
    HYPERBOLIC_CFL,USE_POISSON_SOLVER_5GLL

  use specfem_par, only: SIMULATION_TYPE,deltat,scale_veloc

  use specfem_par_full_gravity, only: neq1

  implicit none

  ! local parameters
  integer :: ier
  double precision :: cg,vp_max,r_outer,lambda_max,dt_stable

  ! the damped equation is not time-reversible: run backwards the damping term
  ! amplifies instead of relaxing, so the reconstructed forward wavefield of a
  ! kernel simulation would blow up
  if (SIMULATION_TYPE == 3) &
    call exit_MPI(myrank,'POISSON_SOLVER = 2 (hyperbolic) does not support SIMULATION_TYPE = 3 yet')

  ! only the Level-1 system is marched
  if (USE_POISSON_SOLVER_5GLL) &
    call exit_MPI(myrank,'POISSON_SOLVER = 2 (hyperbolic) requires USE_POISSON_SOLVER_5GLL = .false.')

  ! the Level-1 points must be a subset of the Level-2 points for the geometry
  ! of the elastic mesh to be reusable by strided access
  if (NGLLX_INF < 2 .or. mod(NGLLX-1,NGLLX_INF-1) /= 0) &
    call exit_MPI(myrank,'POISSON_SOLVER = 2 (hyperbolic) needs NGLLX-1 divisible by NGLLX_INF-1')

  igll_stride = (NGLLX-1) / (NGLLX_INF-1)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    initializing HYPERBOLIC solver'
    call flush_IMAIN()
  endif

  allocate(phi_dot(0:neq1),phi_ddot(0:neq1), &
           minv1(0:neq1),load_prev(0:neq1), &
           work_res(0:neq1),work_acc(0:neq1),stat=ier)
  if (ier /= 0) stop 'Error allocating hyperbolic solver arrays'

  phi_dot(:)   = 0.0_CUSTOM_REAL
  phi_ddot(:)  = 0.0_CUSTOM_REAL
  minv1(:)     = 0.0_CUSTOM_REAL
  load_prev(:) = 0.0_CUSTOM_REAL
  work_res(:)  = 0.0_CUSTOM_REAL
  work_acc(:)  = 0.0_CUSTOM_REAL

  ! artificial gravity speed, referenced to the fastest seismic wave on the mesh so
  ! that no propagating mode is left under-resolved
  call get_max_pvelocity(vp_max)

  cg = HYPERBOLIC_KG * vp_max
  cg2 = real(cg*cg, kind=CUSTOM_REAL)

  ! 1D GLL operators, metric factors and lumped mass
  call build_gll_operators()
  call build_operator_and_mass()

  ! largest eigenvalue of M^-1 K sets the stable sub-step
  call estimate_lambda_max(lambda_max)

  ! explicit Newmark is stable for dt <= 2/omega_max with omega_max = cg*sqrt(lambda_max);
  ! the 1.1 covers the power iteration approaching lambda_max from below
  dt_stable = HYPERBOLIC_CFL * 2.d0 / (cg * dsqrt(1.1d0 * lambda_max))

  nsub = max(1, ceiling(dble(deltat) / dt_stable))
  dt_sub = deltat / real(nsub, kind=CUSTOM_REAL)

  ! damping: kappa = kappa_tilde * cg / L, so that kappa_tilde ~ pi critically
  ! damps the gravest artificial mode of the finite domain
  call get_outer_radius(r_outer)
  kappa = real(HYPERBOLIC_KAPPA * cg / r_outer, kind=CUSTOM_REAL)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '      Level-1 GLL points        = ',NGLLX_INF,' (stride',igll_stride,'in the elastic mesh)'
    write(IMAIN,*) '      max P velocity of mesh    = ',sngl(vp_max*scale_veloc/1000.d0),'km/s'
    write(IMAIN,*) '      gravity speed cg          = ',sngl(cg*scale_veloc/1000.d0),'km/s'
    write(IMAIN,*) '                                = ',sngl(cg),'(non-dim)'
    write(IMAIN,*) '      relative error bound      = ',sngl(1.d0/(HYPERBOLIC_KG*HYPERBOLIC_KG))
    write(IMAIN,*) '      damping rate kappa        = ',sngl(dble(kappa)),'(non-dim)'
    write(IMAIN,*) '      outer radius of domain    = ',sngl(r_outer),'(non-dim)'
    write(IMAIN,*) '      max eigenvalue of M^-1 K  = ',sngl(lambda_max)
    write(IMAIN,*) '      stable gravity sub-step   = ',sngl(dt_stable)
    write(IMAIN,*) '      sub-steps per elastic step= ',nsub
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! the stored element stiffness matrices of the finite regions are now dead
  ! weight: the operator below builds their action on the fly
  call release_stored_stiffness()

  call synchronize_all()

  end subroutine hyperbolic_prepare

!
!-------------------------------------------------------------------------------
!

  subroutine hyperbolic_solve()

! advances the perturbed potential over one elastic time step

  use constants_solver, only: NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE, &
    NGLOB_CRUST_MANTLE,NGLOB_OUTER_CORE,NGLOB_INNER_CORE,NGLOB_TRINFINITE,NGLOB_INFINITE

  use specfem_par_full_gravity, only: gravload1, &
    pgrav1,pgrav_ic1,pgrav_oc1,pgrav_cm1,pgrav_trinf1,pgrav_inf1, &
    pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf, &
    gdof_ic1,gdof_oc1,gdof_cm1,gdof_trinf1,gdof_inf1, &
    inode_elmt_ic,inode_elmt_oc,inode_elmt_cm,inode_elmt_trinf,inode_elmt_inf, &
    inode_map_ic,inode_map_oc,inode_map_cm,inode_map_trinf,inode_map_inf, &
    nmir_ic,nmir_oc,nmir_cm,nmir_trinf,nmir_inf, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    is_active_gll,igll_active_on

  use siem_poisson, only: compute_poisson_load3

  use siem_solver_mpi, only: interpolate3to5

  implicit none

  ! local parameters
  integer :: isub
  real(kind=CUSTOM_REAL) :: theta

  ! load at the end of this elastic step, from the updated displacement
  call compute_poisson_load3()

  ! sub-cycle the gravity equation across the elastic step, ramping the source
  ! linearly from the previous step's load to this one
  do isub = 1,nsub
    theta = real(isub, kind=CUSTOM_REAL) / real(nsub, kind=CUSTOM_REAL)
    call hyperbolic_advance(theta)
  enddo

  load_prev(:) = gravload1(:)

  ! interpolate the Level-1 potential back onto the Level-2 mesh that the
  ! elastic and fluid solvers read
  pgrav_ic1(:) = pgrav1(gdof_ic1(:))
  pgrav_oc1(:) = pgrav1(gdof_oc1(:))
  pgrav_cm1(:) = pgrav1(gdof_cm1(:))
  if (ADD_TRINF) pgrav_trinf1(:) = pgrav1(gdof_trinf1(:))
  pgrav_inf1(:) = pgrav1(gdof_inf1(:))

  call interpolate3to5(NSPEC_INNER_CORE,NGLOB_INNER_CORE,nnode_ic1, &
                       inode_elmt_ic,nmir_ic,inode_map_ic,is_active_gll,igll_active_on,pgrav_ic1,pgrav_ic)

  call interpolate3to5(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,nnode_oc1, &
                       inode_elmt_oc,nmir_oc,inode_map_oc,is_active_gll,igll_active_on,pgrav_oc1,pgrav_oc)

  call interpolate3to5(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,nnode_cm1, &
                       inode_elmt_cm,nmir_cm,inode_map_cm,is_active_gll,igll_active_on,pgrav_cm1,pgrav_cm)

  if (ADD_TRINF) then
    call interpolate3to5(NSPEC_TRINFINITE,NGLOB_TRINFINITE,nnode_trinf1, &
                         inode_elmt_trinf,nmir_trinf,inode_map_trinf,is_active_gll,igll_active_on,pgrav_trinf1,pgrav_trinf)
  endif

  call interpolate3to5(NSPEC_INFINITE,NGLOB_INFINITE,nnode_inf1, &
                       inode_elmt_inf,nmir_inf,inode_map_inf,is_active_gll,igll_active_on,pgrav_inf1,pgrav_inf)

  end subroutine hyperbolic_solve

!
!-------------------------------------------------------------------------------
!

  subroutine hyperbolic_advance(theta)

! one explicit Newmark sub-step of  M phi_tt + 2 kappa M phi_t + cg^2 K phi = cg^2 F
!
! this is the scalar counterpart of update_displ_Newmark() followed by
! compute_forces_*() and the velocity/acceleration correction

  use specfem_par_full_gravity, only: neq1,gravload1,pgrav1

  use siem_solver_mpi, only: scatter_and_assemble3

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: theta

  ! local parameters
  real(kind=CUSTOM_REAL) :: dt_half,dt_sq_half,damp

  dt_half = 0.5_CUSTOM_REAL * dt_sub
  dt_sq_half = dt_half * dt_sub

  ! predictor: potential to the end of the sub-step, velocity to the mid-point
  pgrav1(:) = pgrav1(:) + dt_sub * phi_dot(:) + dt_sq_half * phi_ddot(:)
  phi_dot(:) = phi_dot(:) + dt_half * phi_ddot(:)

  pgrav1(0) = 0.0_CUSTOM_REAL
  phi_dot(0) = 0.0_CUSTOM_REAL

  ! residual cg^2 * (F - K phi), regionally assembled only, matching gravload1
  call gravity_stiffness(neq1,pgrav1,work_res)

  work_res(:) = cg2 * ( ((1.0_CUSTOM_REAL - theta) * load_prev(:) + theta * gravload1(:)) - work_res(:) )

  call scatter_and_assemble3(neq1,work_res,work_acc)

  ! undamped acceleration
  work_acc(:) = minv1(:) * work_acc(:)

  ! corrector. The damping term is linear in the end-of-step velocity, so it is
  ! resolved exactly rather than lagged, which keeps large kappa stable
  damp = 1.0_CUSTOM_REAL / (1.0_CUSTOM_REAL + kappa * dt_sub)

  phi_dot(:) = damp * (phi_dot(:) + dt_half * work_acc(:))
  phi_ddot(:) = work_acc(:) - 2.0_CUSTOM_REAL * kappa * phi_dot(:)

  phi_dot(0) = 0.0_CUSTOM_REAL
  phi_ddot(0) = 0.0_CUSTOM_REAL

  end subroutine hyperbolic_advance

!
!-------------------------------------------------------------------------------
!

  subroutine gravity_stiffness(neq,phi_g,kphi)

! action of the Poisson stiffness on the potential, K phi.
!
! Same contract as the SIEM routine it replaces: phi_g must be fully assembled,
! and kphi comes back assembled across regions but not across MPI processes.
! The finite regions are matrix-free; the infinite layer keeps its stored
! element matrices because its shape functions are not a tensor-product basis.

  use specfem_par, only: NSPEC_INNER_CORE,NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE

  use specfem_par_full_gravity, only: &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    gdof_cm1,inode_elmt_cm1, &
    gdof_oc1,inode_elmt_oc1, &
    gdof_ic1,inode_elmt_ic1, &
    gdof_trinf1,inode_elmt_trinf1, &
    gdof_inf1,inode_elmt_inf1,storekmat_infinite1

  implicit none

  integer,intent(in) :: neq
  real(kind=CUSTOM_REAL),intent(in) :: phi_g(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: kphi(0:neq)

  ! local parameters
  real(kind=CUSTOM_REAL) :: kp_ic(nnode_ic1),kp_oc(nnode_oc1),kp_cm(nnode_cm1), &
                            kp_trinf(nnode_trinf1),kp_inf(nnode_inf1)
  real(kind=CUSTOM_REAL) :: km_inf(NGLLCUBE_INF,NGLLCUBE_INF)
  integer :: i,i_elmt
  integer :: inode_inf(NGLLCUBE_INF),igdof_inf(NGLLCUBE_INF)

  real(kind=CUSTOM_REAL),parameter :: zero = 0.0_CUSTOM_REAL

  ! inner core
  call region_stiffness(NSPEC_INNER_CORE,nnode_ic1,inode_elmt_ic1,gdof_ic1,gfac_ic,neq,phi_g,kp_ic)

  ! outer core
  call region_stiffness(NSPEC_OUTER_CORE,nnode_oc1,inode_elmt_oc1,gdof_oc1,gfac_oc,neq,phi_g,kp_oc)

  ! crust mantle
  call region_stiffness(NSPEC_CRUST_MANTLE,nnode_cm1,inode_elmt_cm1,gdof_cm1,gfac_cm,neq,phi_g,kp_cm)

  ! transition infinite
  kp_trinf(:) = zero
  if (ADD_TRINF) then
    call region_stiffness(NSPEC_TRINFINITE,nnode_trinf1,inode_elmt_trinf1,gdof_trinf1,gfac_trinf, &
                          neq,phi_g,kp_trinf)
  endif

  ! infinite layer: mapped-infinite elements, kept matrix-based
  kp_inf(:) = zero
  do i_elmt = 1,NSPEC_INFINITE
    inode_inf(:) = inode_elmt_inf1(:,i_elmt)
    igdof_inf(:) = gdof_inf1(inode_inf(:))
    km_inf = storekmat_infinite1(:,:,i_elmt)
    kp_inf(inode_inf(:)) = kp_inf(inode_inf(:)) + matmul(km_inf,phi_g(igdof_inf))
  enddo

  ! assemble across the regions of this process, but not across MPI
  kphi(:) = zero

  do i = 1,nnode_cm1
    kphi(gdof_cm1(i)) = kphi(gdof_cm1(i)) + kp_cm(i)
  enddo
  do i = 1,nnode_oc1
    kphi(gdof_oc1(i)) = kphi(gdof_oc1(i)) + kp_oc(i)
  enddo
  do i = 1,nnode_ic1
    kphi(gdof_ic1(i)) = kphi(gdof_ic1(i)) + kp_ic(i)
  enddo
  if (ADD_TRINF) then
    do i = 1,nnode_trinf1
      kphi(gdof_trinf1(i)) = kphi(gdof_trinf1(i)) + kp_trinf(i)
    enddo
  endif
  do i = 1,nnode_inf1
    kphi(gdof_inf1(i)) = kphi(gdof_inf1(i)) + kp_inf(i)
  enddo

  kphi(0) = zero

  end subroutine gravity_stiffness

!
!-------------------------------------------------------------------------------
!

  subroutine region_stiffness(nelmt,nnode1,inode_elmt1,gdof1,gfac,neq,phi_g,kp)

! matrix-free K phi over one region, gathering from and scattering to the
! region's own node numbering

  implicit none

  integer,intent(in) :: nelmt,nnode1,neq
  integer,intent(in) :: inode_elmt1(NGLLCUBE_INF,nelmt)
  integer,intent(in) :: gdof1(nnode1)
  real(kind=CUSTOM_REAL),intent(in) :: gfac(6,NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: phi_g(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: kp(nnode1)

  ! local parameters
  integer :: i,j,k,igll,i_elmt
  integer :: inode(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: u(NGLLX_INF,NGLLY_INF,NGLLZ_INF)
  real(kind=CUSTOM_REAL) :: ku(NGLLX_INF,NGLLY_INF,NGLLZ_INF)

  kp(:) = 0.0_CUSTOM_REAL

  do i_elmt = 1,nelmt
    inode(:) = inode_elmt1(:,i_elmt)

    ! gather. inode_elmt1 is ordered with x fastest, then y, then z
    do k = 1,NGLLZ_INF
      do j = 1,NGLLY_INF
        do i = 1,NGLLX_INF
          igll = NGLLX_INF*NGLLY_INF*(k-1) + NGLLX_INF*(j-1) + i
          u(i,j,k) = phi_g(gdof1(inode(igll)))
        enddo
      enddo
    enddo

    call element_laplacian(gfac(:,:,i_elmt),u,ku)

    do k = 1,NGLLZ_INF
      do j = 1,NGLLY_INF
        do i = 1,NGLLX_INF
          igll = NGLLX_INF*NGLLY_INF*(k-1) + NGLLX_INF*(j-1) + i
          kp(inode(igll)) = kp(inode(igll)) + ku(i,j,k)
        enddo
      enddo
    enddo
  enddo

  end subroutine region_stiffness

!
!-------------------------------------------------------------------------------
!

  subroutine element_laplacian(gfac_e,u,ku)

! element action of the weak Laplacian, K_ij = integral( grad N_i . grad N_j ).
!
! Tensor-product form, identical in structure to the derivative contractions in
! compute_forces_outer_core_Dev: differentiate in each reference direction,
! raise the index with the metric, then contract back against the weighted
! derivative of the test functions.

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: gfac_e(6,NGLLCUBE_INF)
  real(kind=CUSTOM_REAL),intent(in) :: u(NGLLX_INF,NGLLY_INF,NGLLZ_INF)
  real(kind=CUSTOM_REAL),intent(out) :: ku(NGLLX_INF,NGLLY_INF,NGLLZ_INF)

  ! local parameters
  real(kind=CUSTOM_REAL) :: t1(NGLLX_INF,NGLLY_INF,NGLLZ_INF)
  real(kind=CUSTOM_REAL) :: t2(NGLLX_INF,NGLLY_INF,NGLLZ_INF)
  real(kind=CUSTOM_REAL) :: t3(NGLLX_INF,NGLLY_INF,NGLLZ_INF)
  real(kind=CUSTOM_REAL) :: du1,du2,du3,s1,s2,s3
  integer :: i,j,k,l,igll

  ! reference-space gradient, raised to the contravariant flux
  do k = 1,NGLLZ_INF
    do j = 1,NGLLY_INF
      do i = 1,NGLLX_INF
        du1 = 0.0_CUSTOM_REAL
        du2 = 0.0_CUSTOM_REAL
        du3 = 0.0_CUSTOM_REAL
        do l = 1,NGLLX_INF
          du1 = du1 + hprime_gll(i,l) * u(l,j,k)
        enddo
        do l = 1,NGLLY_INF
          du2 = du2 + hprime_gll(j,l) * u(i,l,k)
        enddo
        do l = 1,NGLLZ_INF
          du3 = du3 + hprime_gll(k,l) * u(i,j,l)
        enddo

        igll = NGLLX_INF*NGLLY_INF*(k-1) + NGLLX_INF*(j-1) + i

        t1(i,j,k) = gfac_e(1,igll)*du1 + gfac_e(2,igll)*du2 + gfac_e(3,igll)*du3
        t2(i,j,k) = gfac_e(2,igll)*du1 + gfac_e(4,igll)*du2 + gfac_e(5,igll)*du3
        t3(i,j,k) = gfac_e(3,igll)*du1 + gfac_e(5,igll)*du2 + gfac_e(6,igll)*du3
      enddo
    enddo
  enddo

  ! contract against the weighted test-function derivatives
  do k = 1,NGLLZ_INF
    do j = 1,NGLLY_INF
      do i = 1,NGLLX_INF
        s1 = 0.0_CUSTOM_REAL
        s2 = 0.0_CUSTOM_REAL
        s3 = 0.0_CUSTOM_REAL
        do l = 1,NGLLX_INF
          s1 = s1 + hprimew_gll(l,i) * t1(l,j,k)
        enddo
        do l = 1,NGLLY_INF
          s2 = s2 + hprimew_gll(l,j) * t2(i,l,k)
        enddo
        do l = 1,NGLLZ_INF
          s3 = s3 + hprimew_gll(l,k) * t3(i,j,l)
        enddo

        ku(i,j,k) = wgll1d(j)*wgll1d(k)*s1 &
                  + wgll1d(i)*wgll1d(k)*s2 &
                  + wgll1d(i)*wgll1d(j)*s3
      enddo
    enddo
  enddo

  end subroutine element_laplacian

!
!-------------------------------------------------------------------------------
!

  subroutine build_gll_operators()

! 1D GLL points, weights and derivative matrix of the Level-1 basis,
! the analogue of hprime_xx / hprimewgll_xx for NGLLX_INF points

  use siem_gll_library, only: kdble,zwgljd,lagrange1dGLLAS

  implicit none

  ! local parameters
  real(kind=kdble),parameter :: jalpha = 0.0_kdble,jbeta = 0.0_kdble
  real(kind=kdble) :: xigll(NGLLX_INF),wgll(NGLLX_INF)
  real(kind=kdble) :: phi(NGLLX_INF),dphi(NGLLX_INF)
  integer :: i,l

  call zwgljd(xigll,wgll,NGLLX_INF,jalpha,jbeta)

  do i = 1,NGLLX_INF
    wgll1d(i) = real(wgll(i), kind=CUSTOM_REAL)
  enddo

  do i = 1,NGLLX_INF
    call lagrange1dGLLAS(NGLLX_INF,xigll,xigll(i),phi,dphi)
    do l = 1,NGLLX_INF
      hprime_gll(i,l) = real(dphi(l), kind=CUSTOM_REAL)
    enddo
  enddo

  ! weighted transpose used by the second contraction
  do i = 1,NGLLX_INF
    do l = 1,NGLLX_INF
      hprimew_gll(l,i) = hprime_gll(l,i) * wgll1d(l)
    enddo
  enddo

  end subroutine build_gll_operators

!
!-------------------------------------------------------------------------------
!

  subroutine build_operator_and_mass()

! metric factors of the matrix-free Laplacian and the lumped mass matrix,
! both taken from the geometry the elastic solver already stores

  use constants, only: IREGION_INNER_CORE,IREGION_OUTER_CORE,IREGION_CRUST_MANTLE, &
    IREGION_TRINFINITE,HYPERBOLIC_MASS_INF

  use constants_solver, only: NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE, &
    NSPEC_TRINFINITE

  use specfem_par, only: NPROCTOT_VAL

  use specfem_par_crustmantle, only: xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  use specfem_par_outercore, only: xix_outer_core,xiy_outer_core,xiz_outer_core, &
    etax_outer_core,etay_outer_core,etaz_outer_core, &
    gammax_outer_core,gammay_outer_core,gammaz_outer_core
  use specfem_par_innercore, only: xix_inner_core,xiy_inner_core,xiz_inner_core, &
    etax_inner_core,etay_inner_core,etaz_inner_core, &
    gammax_inner_core,gammay_inner_core,gammaz_inner_core
  use specfem_par_trinfinite, only: xix_trinfinite,xiy_trinfinite,xiz_trinfinite, &
    etax_trinfinite,etay_trinfinite,etaz_trinfinite, &
    gammax_trinfinite,gammay_trinfinite,gammaz_trinfinite

  use specfem_par_full_gravity, only: neq1,dprecon1, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    inode_elmt_ic1,inode_elmt_oc1,inode_elmt_cm1,inode_elmt_trinf1, &
    gdof_ic1,gdof_oc1,gdof_cm1,gdof_trinf1,gdof_inf1, &
    num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
    nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1,my_neighbors_crust_mantle1, &
    num_interfaces_outer_core1,max_nibool_interfaces_outer_core1, &
    nibool_interfaces_outer_core1,ibool_interfaces_outer_core1,my_neighbors_outer_core1, &
    num_interfaces_inner_core1,max_nibool_interfaces_inner_core1, &
    nibool_interfaces_inner_core1,ibool_interfaces_inner_core1,my_neighbors_inner_core1, &
    num_interfaces_trinfinite1,max_nibool_interfaces_trinfinite1, &
    nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1,my_neighbors_trinfinite1

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: mass_ic,mass_oc,mass_cm,mass_trinf
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: mass1
  logical,dimension(:),allocatable :: is_finite_dof
  real(kind=CUSTOM_REAL) :: mass_sum,stiff_sum,mass_sum_all,stiff_sum_all,ratio
  double precision :: sizeval
  integer :: i,ier

  allocate(gfac_cm(6,NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           gfac_oc(6,NGLLCUBE_INF,NSPEC_OUTER_CORE), &
           gfac_ic(6,NGLLCUBE_INF,NSPEC_INNER_CORE), &
           gfac_trinf(6,NGLLCUBE_INF,NSPEC_TRINFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating hyperbolic metric factor arrays'

  gfac_cm(:,:,:) = 0.0_CUSTOM_REAL
  gfac_oc(:,:,:) = 0.0_CUSTOM_REAL
  gfac_ic(:,:,:) = 0.0_CUSTOM_REAL
  gfac_trinf(:,:,:) = 0.0_CUSTOM_REAL

  allocate(mass_ic(nnode_ic1),mass_oc(nnode_oc1),mass_cm(nnode_cm1), &
           mass_trinf(nnode_trinf1),mass1(0:neq1), &
           is_finite_dof(0:neq1),stat=ier)
  if (ier /= 0) stop 'Error allocating hyperbolic mass arrays'

  call region_geometry(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,nnode_cm1,inode_elmt_cm1, &
                       xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                       etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                       gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                       gfac_cm,mass_cm)

  call region_geometry(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,nnode_oc1,inode_elmt_oc1, &
                       xix_outer_core,xiy_outer_core,xiz_outer_core, &
                       etax_outer_core,etay_outer_core,etaz_outer_core, &
                       gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                       gfac_oc,mass_oc)

  call region_geometry(IREGION_INNER_CORE,NSPEC_INNER_CORE,nnode_ic1,inode_elmt_ic1, &
                       xix_inner_core,xiy_inner_core,xiz_inner_core, &
                       etax_inner_core,etay_inner_core,etaz_inner_core, &
                       gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                       gfac_ic,mass_ic)

  mass_trinf(:) = 0.0_CUSTOM_REAL
  if (ADD_TRINF) then
    call region_geometry(IREGION_TRINFINITE,NSPEC_TRINFINITE,nnode_trinf1,inode_elmt_trinf1, &
                         xix_trinfinite,xiy_trinfinite,xiz_trinfinite, &
                         etax_trinfinite,etay_trinfinite,etaz_trinfinite, &
                         gammax_trinfinite,gammay_trinfinite,gammaz_trinfinite, &
                         gfac_trinf,mass_trinf)
  endif

  ! assemble the mass across MPI processes within each region, then across
  ! regions, following the same protocol as the stiffness preconditioner
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_cm1,mass_cm, &
                           num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
                           nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1, &
                           my_neighbors_crust_mantle1)

  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_oc1,mass_oc, &
                           num_interfaces_outer_core1,max_nibool_interfaces_outer_core1, &
                           nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
                           my_neighbors_outer_core1)

  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_ic1,mass_ic, &
                           num_interfaces_inner_core1,max_nibool_interfaces_inner_core1, &
                           nibool_interfaces_inner_core1,ibool_interfaces_inner_core1, &
                           my_neighbors_inner_core1)

  if (ADD_TRINF) then
    call assemble_MPI_scalar(NPROCTOT_VAL,nnode_trinf1,mass_trinf, &
                             num_interfaces_trinfinite1,max_nibool_interfaces_trinfinite1, &
                             nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1, &
                             my_neighbors_trinfinite1)
  endif

  mass1(:) = 0.0_CUSTOM_REAL
  is_finite_dof(:) = .false.

  do i = 1,nnode_cm1
    mass1(gdof_cm1(i)) = mass1(gdof_cm1(i)) + mass_cm(i)
    is_finite_dof(gdof_cm1(i)) = .true.
  enddo
  do i = 1,nnode_oc1
    mass1(gdof_oc1(i)) = mass1(gdof_oc1(i)) + mass_oc(i)
    is_finite_dof(gdof_oc1(i)) = .true.
  enddo
  do i = 1,nnode_ic1
    mass1(gdof_ic1(i)) = mass1(gdof_ic1(i)) + mass_ic(i)
    is_finite_dof(gdof_ic1(i)) = .true.
  enddo
  if (ADD_TRINF) then
    do i = 1,nnode_trinf1
      mass1(gdof_trinf1(i)) = mass1(gdof_trinf1(i)) + mass_trinf(i)
      is_finite_dof(gdof_trinf1(i)) = .true.
    enddo
  endif

  is_finite_dof(0) = .false.

  ! the infinite elements have no finite mass integral, so their interior
  ! degrees of freedom get a pseudo-mass carrying the same mass-to-stiffness
  ! ratio as the outermost finite region. Those elements are the largest in the
  ! mesh, so this keeps the exterior off the critical path of the CFL while
  ! letting it relax to the Laplace solution within a few sub-steps. The fixed
  ! point K phi = F does not depend on M, so this cannot bias the answer.
  if (ADD_TRINF) then
    mass_sum = sum(mass_trinf)
    stiff_sum = sum(dprecon1(gdof_trinf1(:)))
  else
    mass_sum = sum(mass_cm)
    stiff_sum = sum(dprecon1(gdof_cm1(:)))
  endif

  call sum_all_all_cr(mass_sum,mass_sum_all)
  call sum_all_all_cr(stiff_sum,stiff_sum_all)

  if (stiff_sum_all <= 0.0_CUSTOM_REAL) &
    call exit_MPI(myrank,'Invalid stiffness sum while building hyperbolic pseudo-mass')

  ratio = real(HYPERBOLIC_MASS_INF, kind=CUSTOM_REAL) * mass_sum_all / stiff_sum_all

  do i = 1,nnode_inf1
    if (gdof_inf1(i) > 0) then
      if (.not. is_finite_dof(gdof_inf1(i))) then
        mass1(gdof_inf1(i)) = ratio * dprecon1(gdof_inf1(i))
      endif
    endif
  enddo

  ! invert
  minv1(:) = 0.0_CUSTOM_REAL
  do i = 1,neq1
    if (mass1(i) <= 0.0_CUSTOM_REAL) &
      call exit_MPI(myrank,'Zero mass on an active degree of freedom in the hyperbolic solver')
    minv1(i) = 1.0_CUSTOM_REAL / mass1(i)
  enddo

  ! user output
  if (myrank == 0) then
    sizeval = 6.d0 * dble(NGLLCUBE_INF) * dble(CUSTOM_REAL) &
              * (dble(NSPEC_CRUST_MANTLE) + dble(NSPEC_OUTER_CORE) + dble(NSPEC_INNER_CORE) &
                 + dble(NSPEC_TRINFINITE)) / 1024.d0 / 1024.d0
    write(IMAIN,*) '      matrix-free metric factors= ',sngl(sizeval),'MB'
    call flush_IMAIN()
  endif

  deallocate(mass_ic,mass_oc,mass_cm,mass_trinf,mass1,is_finite_dof)

  end subroutine build_operator_and_mass

!
!-------------------------------------------------------------------------------
!

  subroutine region_geometry(iregion,nelmt,nnode1,inode_elmt1, &
                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                             gfac,mass)

! metric factors and lumped mass of one region, read off the elastic mesh
! geometry at the Level-1 subset of the GLL points

  use constants_solver, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE

  use specfem_par_innercore, only: idoubling_inner_core

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode1
  integer,intent(in) :: inode_elmt1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nelmt),intent(in) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL),intent(out) :: gfac(6,NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(out) :: mass(nnode1)

  ! local parameters
  integer :: i,j,k,i5,j5,k5,igll,i_elmt
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: detinv,jacobianl

  gfac(:,:,:) = 0.0_CUSTOM_REAL
  mass(:) = 0.0_CUSTOM_REAL

  do i_elmt = 1,nelmt
    ! fictitious elements of the central cube carry no stiffness or mass
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    endif

    do k = 1,NGLLZ_INF
      k5 = (k-1)*igll_stride + 1
      do j = 1,NGLLY_INF
        j5 = (j-1)*igll_stride + 1
        do i = 1,NGLLX_INF
          i5 = (i-1)*igll_stride + 1

          xixl = xix(i5,j5,k5,i_elmt)
          xiyl = xiy(i5,j5,k5,i_elmt)
          xizl = xiz(i5,j5,k5,i_elmt)
          etaxl = etax(i5,j5,k5,i_elmt)
          etayl = etay(i5,j5,k5,i_elmt)
          etazl = etaz(i5,j5,k5,i_elmt)
          gammaxl = gammax(i5,j5,k5,i_elmt)
          gammayl = gammay(i5,j5,k5,i_elmt)
          gammazl = gammaz(i5,j5,k5,i_elmt)

          ! stored arrays hold the inverse mapping, so this determinant inverts
          detinv = xixl*(etayl*gammazl - etazl*gammayl) &
                 - xiyl*(etaxl*gammazl - etazl*gammaxl) &
                 + xizl*(etaxl*gammayl - etayl*gammaxl)

          if (detinv == 0.0_CUSTOM_REAL) &
            call exit_MPI(myrank,'Singular Jacobian while building the hyperbolic operator')

          jacobianl = 1.0_CUSTOM_REAL / detinv

          igll = NGLLX_INF*NGLLY_INF*(k-1) + NGLLX_INF*(j-1) + i

          ! J * contravariant metric, packed as (g11,g12,g13,g22,g23,g33)
          gfac(1,igll,i_elmt) = jacobianl * (xixl*xixl + xiyl*xiyl + xizl*xizl)
          gfac(2,igll,i_elmt) = jacobianl * (xixl*etaxl + xiyl*etayl + xizl*etazl)
          gfac(3,igll,i_elmt) = jacobianl * (xixl*gammaxl + xiyl*gammayl + xizl*gammazl)
          gfac(4,igll,i_elmt) = jacobianl * (etaxl*etaxl + etayl*etayl + etazl*etazl)
          gfac(5,igll,i_elmt) = jacobianl * (etaxl*gammaxl + etayl*gammayl + etazl*gammazl)
          gfac(6,igll,i_elmt) = jacobianl * (gammaxl*gammaxl + gammayl*gammayl + gammazl*gammazl)

          ! GLL collocation makes the element mass matrix diagonal
          mass(inode_elmt1(igll,i_elmt)) = mass(inode_elmt1(igll,i_elmt)) &
                                           + jacobianl * wgll1d(i) * wgll1d(j) * wgll1d(k)
        enddo
      enddo
    enddo
  enddo

  end subroutine region_geometry

!
!-------------------------------------------------------------------------------
!

  subroutine release_stored_stiffness()

! the finite regions no longer need their stored element stiffness matrices.
! The infinite layer keeps its own, and dprecon1 has already been formed.

  use specfem_par_full_gravity, only: storekmat_crust_mantle1,storekmat_outer_core1, &
    storekmat_inner_core1,storekmat_trinfinite1

  implicit none

  ! local parameters
  double precision :: sizeval

  sizeval = 0.d0

  if (allocated(storekmat_crust_mantle1)) then
    sizeval = sizeval + dble(size(storekmat_crust_mantle1))
    deallocate(storekmat_crust_mantle1)
  endif
  if (allocated(storekmat_outer_core1)) then
    sizeval = sizeval + dble(size(storekmat_outer_core1))
    deallocate(storekmat_outer_core1)
  endif
  if (allocated(storekmat_inner_core1)) then
    sizeval = sizeval + dble(size(storekmat_inner_core1))
    deallocate(storekmat_inner_core1)
  endif
  if (allocated(storekmat_trinfinite1)) then
    sizeval = sizeval + dble(size(storekmat_trinfinite1))
    deallocate(storekmat_trinfinite1)
  endif

  sizeval = sizeval * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0

  if (myrank == 0) then
    write(IMAIN,*) '      released element stiffness= ',sngl(sizeval),'MB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine release_stored_stiffness

!
!-------------------------------------------------------------------------------
!

  subroutine estimate_lambda_max(lambda_max)

! power iteration for the largest eigenvalue of M^-1 K, which sets the CFL

  use constants, only: HYPERBOLIC_NPOWER

  use specfem_par_full_gravity, only: neq1

  use siem_math_library_mpi, only: maxvec

  use siem_solver_mpi, only: scatter_and_assemble3

  implicit none

  double precision,intent(out) :: lambda_max

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: v,kv,av
  real(kind=CUSTOM_REAL) :: vnorm,anorm
  integer :: iter,i,ier

  allocate(v(0:neq1),kv(0:neq1),av(0:neq1),stat=ier)
  if (ier /= 0) stop 'Error allocating hyperbolic power iteration arrays'

  ! a deterministic but non-smooth start, so the iterate is not orthogonal to
  ! the highest mode (which is what a constant vector would be)
  kv(0) = 0.0_CUSTOM_REAL
  do i = 1,neq1
    kv(i) = real(1 - 2*mod(i,2), kind=CUSTOM_REAL)
  enddo

  ! the local ordering differs between ranks, so make the start vector agree
  ! across processes before iterating with it
  call scatter_and_assemble3(neq1,kv,v)

  lambda_max = 0.d0

  do iter = 1,HYPERBOLIC_NPOWER
    call gravity_stiffness(neq1,v,kv)
    call scatter_and_assemble3(neq1,kv,av)

    av(:) = minv1(:) * av(:)
    av(0) = 0.0_CUSTOM_REAL

    vnorm = maxvec(abs(v))
    anorm = maxvec(abs(av))

    if (vnorm <= 0.0_CUSTOM_REAL .or. anorm <= 0.0_CUSTOM_REAL) exit

    lambda_max = dble(anorm) / dble(vnorm)

    v(:) = av(:) / anorm
  enddo

  deallocate(v,kv,av)

  if (lambda_max <= 0.d0) &
    call exit_MPI(myrank,'Power iteration failed to bound M^-1 K in the hyperbolic solver')

  end subroutine estimate_lambda_max

!
!-------------------------------------------------------------------------------
!

  subroutine get_max_pvelocity(vp_max)

! fastest P velocity anywhere on the mesh, in non-dimensional units.
! The isotropic reference moduli are used even for anisotropic models: they are
! always read in, and cg only needs this to within a few percent.

  use specfem_par_crustmantle, only: rhostore_crust_mantle, &
    kappavstore_crust_mantle,muvstore_crust_mantle
  use specfem_par_outercore, only: rhostore_outer_core,kappavstore_outer_core
  use specfem_par_innercore, only: rhostore_inner_core, &
    kappavstore_inner_core,muvstore_inner_core

  implicit none

  double precision,intent(out) :: vp_max

  ! local parameters
  real(kind=CUSTOM_REAL) :: vp2,vp2_all
  real(kind=CUSTOM_REAL),parameter :: FOUR_THIRDS = 4.0_CUSTOM_REAL / 3.0_CUSTOM_REAL

  vp2 = 0.0_CUSTOM_REAL

  ! solid regions: vp^2 = (kappa + 4/3 mu) / rho.
  ! merge() keeps the divisor away from zero where rho is unset (fictitious cube),
  ! since the whole expression is evaluated before the mask selects from it
  vp2 = max(vp2, maxval((kappavstore_crust_mantle + FOUR_THIRDS*muvstore_crust_mantle) &
                        / merge(rhostore_crust_mantle,1.0_CUSTOM_REAL, &
                                rhostore_crust_mantle > 0.0_CUSTOM_REAL), &
                        mask = rhostore_crust_mantle > 0.0_CUSTOM_REAL))

  vp2 = max(vp2, maxval((kappavstore_inner_core + FOUR_THIRDS*muvstore_inner_core) &
                        / merge(rhostore_inner_core,1.0_CUSTOM_REAL, &
                                rhostore_inner_core > 0.0_CUSTOM_REAL), &
                        mask = rhostore_inner_core > 0.0_CUSTOM_REAL))

  ! fluid outer core: vp^2 = kappa / rho
  vp2 = max(vp2, maxval(kappavstore_outer_core &
                        / merge(rhostore_outer_core,1.0_CUSTOM_REAL, &
                                rhostore_outer_core > 0.0_CUSTOM_REAL), &
                        mask = rhostore_outer_core > 0.0_CUSTOM_REAL))

  call max_all_all_cr(vp2,vp2_all)

  if (vp2_all <= 0.0_CUSTOM_REAL) &
    call exit_MPI(myrank,'Could not determine a maximum P velocity for the hyperbolic solver')

  vp_max = dsqrt(dble(vp2_all))

  end subroutine get_max_pvelocity

!
!-------------------------------------------------------------------------------
!

  subroutine get_outer_radius(r_outer)

! outer radius of the finite part of the gravity domain, used to set the
! damping rate relative to the gravest artificial mode it can support

  use constants_solver, only: NGLOB_CRUST_MANTLE,NGLOB_TRINFINITE

  use specfem_par_crustmantle, only: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  use specfem_par_trinfinite, only: xstore_trinfinite,ystore_trinfinite,zstore_trinfinite

  implicit none

  double precision,intent(out) :: r_outer

  ! local parameters
  real(kind=CUSTOM_REAL) :: rmax,rmax_all
  integer :: i

  rmax = 0.0_CUSTOM_REAL

  if (ADD_TRINF) then
    do i = 1,NGLOB_TRINFINITE
      rmax = max(rmax, sqrt(xstore_trinfinite(i)**2 + ystore_trinfinite(i)**2 + zstore_trinfinite(i)**2))
    enddo
  else
    do i = 1,NGLOB_CRUST_MANTLE
      rmax = max(rmax, sqrt(xstore_crust_mantle(i)**2 + ystore_crust_mantle(i)**2 + zstore_crust_mantle(i)**2))
    enddo
  endif

  call max_all_all_cr(rmax,rmax_all)

  r_outer = dble(rmax_all)

  if (r_outer <= 0.d0) &
    call exit_MPI(myrank,'Invalid outer radius for the hyperbolic gravity domain')

  end subroutine get_outer_radius

end module siem_hyperbolic
