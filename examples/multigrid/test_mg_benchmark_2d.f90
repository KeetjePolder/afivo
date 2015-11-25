!> \example test_mg_benchmark_2d.f90
!> Benchmark of the multigrid routines
program test_mg2_2d
  use m_afivo_2d
  use m_mg_2d

  implicit none

  integer, parameter :: dp           = kind(0.0d0)
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 100
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_err = 4

  ! The manufactured solution exists of two Gaussians here.
  ! For each Gaussian, 4 constants are used: pre-factor, x0, y0, sigma.
  integer, parameter :: n_gaussians = 1
  real(dp), parameter :: g_params(4, n_gaussians) = reshape(&
       [1.0_dp, 0.5_dp, 0.5_dp, 0.04_dp], [4, n_gaussians])

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  integer            :: box_size, mesh_size, max_ref_lvl
  real(dp)           :: dr
  character(len=40)  :: var_names(4), arg_string
  type(mg2_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_err) = "err"

  dr = 1.0_dp / box_size

  ! Get box size and mesh size from command line argument
  if (command_argument_count() /= 2) stop "Arguments should be: box_size mesh_size"
  call get_command_argument(1, arg_string)
  read(arg_string, *) box_size
  call get_command_argument(2, arg_string)
  read(arg_string, *) mesh_size

  ! Determine maximum refinement level
  max_ref_lvl = nint(log(mesh_size / real(box_size, dp)) / log(2.0_dp)) + 1

  print *, "Box size: ", box_size
  print *, "Mesh size:", 2**(max_ref_lvl-1) * box_size

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=4, n_var_face=0, &
       dr = dr, coarsen_to = 2, n_boxes=10*1000)

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do
     call a2_loop_box(tree, set_init_cond)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  ! Set the multigrid options
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%n_cycle_down = 2
  mg%n_cycle_up   = 2
  mg%n_cycle_base = 2
  mg%sides_bc     => sides_bc

  call mg2_init_mg(mg)

  do i = 1, n_iterations
     call mg2_fas_fmg(tree, mg, .true., i == 1)
  end do

  call a2_write_silo(tree, "test_mg_benchmark_2d", var_names, 1, 0.0_dp)
  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)

    if (boxes(id)%lvl < max_ref_lvl) ref_flags(id) = a5_do_ref
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_rhs) = rhs(xy)
       end do
    end do
  end subroutine set_init_cond

  subroutine set_err(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - phi_sol(xy)
       end do
    end do
  end subroutine set_err

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    real(dp)                    :: xy(2)
    integer                     :: n, nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [0.5_dp, real(n, dp)])
          boxes(id)%cc(0, n, iv) = 2 * phi_sol(xy) - boxes(id)%cc(1, n, iv)
       end do
    case (a2_nb_hx)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [nc+0.5_dp, real(n, dp)])
          boxes(id)%cc(nc+1, n, iv) = 2 * phi_sol(xy) - boxes(id)%cc(nc, n, iv)
       end do
    case (a2_nb_ly)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [real(n, dp), 0.5_dp])
          boxes(id)%cc(n, 0, iv) = 2 * phi_sol(xy) - boxes(id)%cc(n, 1, iv)
       end do
    case (a2_nb_hy)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [real(n, dp), nc+0.5_dp])
          boxes(id)%cc(n, nc+1, iv) = 2 * phi_sol(xy) - boxes(id)%cc(n, nc, iv)
       end do
    end select
  end subroutine sides_bc

  real(dp) function phi_sol(x)
    real(dp), intent(in) :: x(2)
    integer :: n

    phi_sol = 0
    do n = 1, n_gaussians
       phi_sol = phi_sol + g_params(1, n) * &
            gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function phi_sol

  real(dp) function rhs(x)
    real(dp), intent(in) :: x(2)
    integer :: n

    rhs = 0
    do n = 1, n_gaussians
       rhs = rhs + g_params(1, n) * &
            lpl_gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function rhs

  real(dp) function gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp) :: xrel(2)
    xrel = (x-x0)/sigma
    gaussian_2d = exp(-sum(xrel**2))
  end function gaussian_2d

  real(dp) function lpl_gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp) :: xrel(2)
    xrel = (x-x0)/sigma
    lpl_gaussian_2d = 4/sigma**2 * (sum(xrel**2) - 1) * &
         gaussian_2d(x, x0, sigma)
  end function lpl_gaussian_2d

end program test_mg2_2d
