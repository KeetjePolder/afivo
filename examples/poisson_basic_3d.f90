!> \example poisson_basic_3d.f90
!>
!> Example showing how to use multigrid and compare with an analytic solution. A
!> standard 7-point Laplacian is used.
program poisson_basic_3d
  use m_a3_t
  use m_a3_core
  use m_a3_mg
  use m_a3_utils
  use m_a3_io
  use m_gaussians

  implicit none

  integer, parameter :: n_cell = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 4
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_err = 3
  integer, parameter :: i_tmp = 4

  type(a3_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i
  integer            :: ix_list(3, n_boxes_base)
  integer            :: nb_list(6, n_boxes_base)
  real(dp)           :: dr, min_res, max_res
  character(len=40)  :: fname
  type(mg3_t)        :: mg
  type(gauss_t)      :: gs

  ! The manufactured solution exists of two Gaussians, which are stored in gs
  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([0.25_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.75_dp, 0.5_dp], [3,2]))

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / n_cell

  ! Initialize tree
  call a3_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "err", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1,1]       ! Index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a3_set_base(tree, ix_list, nb_list)

  do
     ! For each box, set the initial conditions
     call a3_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a3_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  call mg3_init_mg(mg)

  do i = 1, 10
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg3_fas_fmg(tree, mg, set_residual=.true., have_guess=(i>1))

     ! Compute the error compared to the analytic solution
     call a3_loop_box(tree, set_error)

     ! Determine the minimum and maximum residual
     call a3_tree_min_cc(tree, i_tmp, min_res)
     call a3_tree_max_cc(tree, i_tmp, max_res)
     print *, "Iteration ", i, "max residual: ", max(abs(min_res), abs(max_res))

     ! This writes a Silo output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "poisson_basic_3d_", i
     call a3_write_silo(tree, trim(fname), dir="output")
  end do

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a3_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    integer                  :: i, j, k, nc
    real(dp)                 :: xyz(3), dr2, drhs

    nc = boxes(id)%n_cell
    dr2 = boxes(id)%dr**2

    outer: do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(boxes(id), [i, j, k])

             ! This is an estimate of the truncation error in the right-hand side,
             ! which is related to the fourth derivative of the solution.
             drhs = dr2 * gauss_4th(gs, xyz) / 12

             if (abs(drhs) > 1.0_dp) then
                ref_flag = a5_do_ref
                exit outer
             end if
          end do
       end do
    end do outer
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             ! Get the coordinate of the cell center at i,j,k
             xyz = a3_r_cc(box, [i,j,k])

             ! And set the rhs values
             box%cc(i, j, k, i_rhs) = gauss_lpl(gs, xyz)
          end do
       end do
    end do
  end subroutine set_init_cond

  ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_err) = box%cc(i, j, k, i_phi) - &
                  gauss_val(gs, xyz)
          end do
       end do
    end do
  end subroutine set_error

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box3_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction in which to set the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    real(dp)                    :: xyz(3), loc
    integer                     :: ix, i, j, k, nc

    nc = box%n_cell

    ! We use dirichlet boundary conditions
    bc_type = a5_bc_dirichlet

    ! Determine whether the direction nb is to "lower" or "higher" neighbors
    if (a3_neighb_low(nb)) then
       ix = 0
       loc = 0.5_dp
    else
       ix = nc+1
       loc = nc + 0.5_dp
    end if

    ! Below the solution is specified in the approriate ghost cells
    select case (a3_neighb_dim(nb))
    case (1)
       do k = 1, nc
          do j = 1, nc
             xyz = a3_rr_cc(box, [loc, real(j, dp), real(k, dp)])
             box%cc(ix, j, k, iv) = gauss_val(gs, xyz)
          end do
       end do
    case (2)
       do k = 1, nc
          do i = 1, nc
             xyz = a3_rr_cc(box, [real(i, dp), loc, real(k, dp)])
             box%cc(i, ix, k, iv) = gauss_val(gs, xyz)
          end do
       end do
    case (3)
       do j = 1, nc
          do i = 1, nc
             xyz = a3_rr_cc(box, [real(i, dp), real(j, dp), loc])
             box%cc(i, j, ix, iv) = gauss_val(gs, xyz)
          end do
       end do
    end select
  end subroutine sides_bc

end program poisson_basic_3d
