!> \example test_mg_cyl_diel.f90
!>
!> Example showing how to use m_a2_mg in cylindrical coordinates with an abrubt
!> change in "eps", and compare with an analytic solution.
program poisson_cyl_dielectric
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_io
  use m_gaussians

  implicit none

  integer, parameter :: n_cell = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 5
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_err = 3
  integer, parameter :: i_tmp = 4
  integer, parameter :: i_eps = 5

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr, max_res, min_res
  character(len=40)  :: fname
  type(mg2_t)        :: mg
  type(gauss_t)      :: gs

  ! The manufactured solution exists of two Gaussians, which are stored in gs
  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp], [2,2]))

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / n_cell

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       coord=a5_cyl, & ! Cylindrical coordinates
       cc_names=["phi", "rhs", "err", "tmp", "eps"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%i_eps        = i_eps       ! Variable for epsilon coefficient
  mg%sides_bc     => sides_bc   ! Method for boundary conditions Because we use

  ! Automatically detect the right methods
  mg%box_op       => mg2_auto_op
  mg%box_gsrb     => mg2_auto_gsrb
  mg%box_corr     => mg2_auto_corr

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary
  call mg2_init_mg(mg)

  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a2_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do

  do i = 1, 10
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, .true., i>1)

     ! Compute the error compared to the analytic solution
     call a2_loop_box(tree, set_err)

     ! Determine the minimum and maximum residual
     call a2_tree_min_cc(tree, i_tmp, min_res)
     call a2_tree_max_cc(tree, i_tmp, max_res)
     print *, "Iteration ", i, "max residual: ", max(abs(min_res), abs(max_res))

     write(fname, "(A,I0)") "poisson_cyl_dielectric_", i
     call a2_write_vtk(tree, trim(fname), dir="output")
  end do

  call a2_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    integer                  :: nc
    real(dp)                 :: max_crv

    nc = boxes(id)%n_cell

    ! Compute the "curvature" in phi
    max_crv = boxes(id)%dr**2 * &
         maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs) / &
         boxes(id)%cc(1:nc, 1:nc, i_eps)))

    ! And refine if it exceeds a threshold
    if (max_crv > 5.0e-4_dp) then
       ref_flag = a5_do_ref
    else
       ref_flag = a5_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2), grad(2), dr, qbnd, tmp

    nc                  = box%n_cell
    box%cc(:, :, i_phi) = 0
    dr                  = box%dr

    do j = 0, nc+1
       do i = 0, nc+1
          rz = a2_r_cc(box, [i,j])

          ! Change epsilon in part of the domain
          if (rz(1) < 0.5_dp .and. rz(2) < 0.5_dp) then
             box%cc(i, j, i_eps) = 100.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if

          ! Partially compute the right-hand side (see below)
          box%cc(i, j, i_rhs) = gauss_lpl_cyl(gs, rz) * box%cc(i, j, i_eps)
       end do
    end do

    ! We have to place surface charges where epsilon has a jump, this is first
    ! done in the r-direction
    do j = 1, nc
       do i = 0, nc
          rz = a2_rr_cc(box, [i + 0.5_dp, real(j, dp)])

          ! Determine amount of charge
          call gauss_grad(gs, rz, grad)
          qbnd = (box%cc(i+1, j, i_eps) - box%cc(i, j, i_eps)) * &
               grad(1) / dr

          ! Place surface charge weighted with eps
          tmp = box%cc(i+1, j, i_eps) / &
               (box%cc(i, j, i_eps) + box%cc(i+1, j, i_eps))
          box%cc(i+1, j, i_rhs) = box%cc(i+1, j, i_rhs) + tmp * qbnd
          box%cc(i, j, i_rhs) = box%cc(i, j, i_rhs) + (1-tmp) * qbnd
       end do
    end do

    ! Set surface charge in z-direction
    do j = 0, nc
       do i = 1, nc
          rz = a2_rr_cc(box, [real(i, dp), j + 0.5_dp])

          ! Determine amount of charge
          call gauss_grad(gs, rz, grad)
          qbnd = (box%cc(i, j+1, i_eps) - box%cc(i, j, i_eps)) * &
               grad(2) / dr

          ! Place surface charge weighted with eps
          tmp = box%cc(i, j+1, i_eps) / &
               (box%cc(i, j, i_eps) + box%cc(i, j+1, i_eps))
          box%cc(i, j+1, i_rhs) = box%cc(i, j+1, i_rhs) + tmp * qbnd
          box%cc(i, j, i_rhs) = box%cc(i, j, i_rhs) + (1-tmp) * qbnd
       end do
    end do
  end subroutine set_init_cond

  ! Compute error with solution
  subroutine set_err(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          rz = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - &
               gauss_val(gs, rz)
       end do
    end do
  end subroutine set_err

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values. Note that on the axis (a boundary in the lower-x
  ! direction) we should use a Neumann zero condition in cylindrical
  ! coordinates.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    real(dp)                    :: rz(2)
    integer                     :: n, nc

    nc = box%n_cell

    select case (nb)
    case (a2_neighb_lowx)             ! Neumann zero on axis
       bc_type = a5_bc_neumann
       box%cc(0, 1:nc, iv) = 0
    case (a2_neighb_highx)             ! Use solution on other boundaries
       bc_type = a5_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = gauss_val(gs, rz)
       end do
    case (a2_neighb_lowy)
       bc_type = a5_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = gauss_val(gs, rz)
       end do
    case (a2_neighb_highy)
       bc_type = a5_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = gauss_val(gs, rz)
       end do
    end select
  end subroutine sides_bc

end program poisson_cyl_dielectric
