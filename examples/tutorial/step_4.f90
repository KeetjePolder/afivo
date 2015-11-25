
program tutorial
  use m_afivo_2d

  implicit none
  integer, parameter  :: dp           = kind(0.0d0)
  integer, parameter  :: box_size     = 8
  integer, parameter  :: n_var_cell   = 1
  integer, parameter  :: n_var_face   = 1
  integer, parameter  :: i_phi        = 1

  real(dp), parameter :: domain_len   = 1.0_dp
  real(dp), parameter :: diff_coeff   = 1.0_dp

  type(a2_t)          :: tree
  integer             :: output_cnt
  real(dp)            :: dt, dt_output, time, end_time
  character(len=100)  :: fname
  type(ref_info_t)    :: ref_info

  print *, "This is an Afivo tutorial for solving the heat equation."

  call initialize_tree(tree)

  do
     call a2_loop_box(tree, set_init_cond)
     call a2_gc_sides(tree, i_phi, a2_sides_interp_lim, a2_bc_dirichlet)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  output_cnt = 0
  time       = 0.0_dp
  end_time   = 1.0_dp
  dt_output  = 0.02_dp

  do while (time < end_time)
     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "heat_", output_cnt
        call a2_write_vtk(tree, fname, ["phi"], output_cnt, time)
     end if

     dt = a2_min_dr(tree)**2 / (6 * diff_coeff)

     call heat_forward_euler(tree, dt)
     time = time + dt

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     call prolong_to_new_children(tree, ref_info)
     call a2_gc_sides(tree, i_phi, a2_sides_interp, a2_bc_dirichlet)
     call a2_tidy_up(tree, 0.8_dp, 0.5_dp, 10000, .false.)
  end do

contains

  subroutine initialize_tree(tree)
    type(a2_t), intent(inout) :: tree

    integer, parameter  :: n_boxes_base = 3
    real(dp)            :: dr
    integer             :: ix_list(2, n_boxes_base)
    integer             :: nb_list(4, n_boxes_base)

    dr = domain_len / box_size

    ! Initialize the tree
    call a2_init(tree, & ! Tree to initialize
         box_size, &     ! Number of cells per coordinate in a box
         n_var_cell, &   ! Number of face-centered variables
         n_var_face, &   ! Number of cell-centered variables
         dr)             ! Distance between cells on base level

    ! Set the indices for boxes 1 to 3
    ix_list(:, 1) = [1, 1]
    ix_list(:, 2) = [2, 1]
    ix_list(:, 3) = [2, 2]

    ! Negative values are used for boundary conditions
    nb_list(1:4, 1) = -1
    nb_list(1:4, 2) = -1
    nb_list(1:4, 3) = -1

    ! Connect the two boxes, nb_hx stands for "neighbor in high-x direction".
    ! Connections only have to be specified from one side.
    nb_list(a2_nb_hx, 1) = 2
    nb_list(a2_nb_hy, 2) = 3

    ! Construct the base (level 1) mesh
    call a2_set_base(tree, ix_list, nb_list)

    call a2_print_info(tree)
  end subroutine initialize_tree

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy-0.5_dp) < 0.25) then
             box%cc(i, j, i_phi) = 1
          else
             box%cc(i, j, i_phi) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine heat_forward_euler(tree, dt)
    type(a2_t), intent(inout) :: tree
    real(dp), intent(in) :: dt

    ! Forward Euler
    call a2_loop_boxes(tree, fluxes_centdif, .true.)
    call a2_consistent_fluxes(tree, [i_phi])
    call a2_loop_box_arg(tree, update_solution, [dt], .true.)
    call a2_restrict_tree(tree, i_phi)
    call a2_gc_sides(tree, i_phi, a2_sides_interp, a2_bc_dirichlet)
  end subroutine heat_forward_euler

  subroutine fluxes_centdif(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    boxes(id)%fx(:,:,i_phi) = (boxes(id)%cc(0:nc, 1:nc, i_phi) - &
         boxes(id)%cc(1:nc+1, 1:nc, i_phi)) * diff_coeff * inv_dr
    boxes(id)%fy(:,:,i_phi) = (boxes(id)%cc(1:nc, 0:nc, i_phi) - &
         boxes(id)%cc(1:nc, 1:nc+1, i_phi)) * diff_coeff * inv_dr
  end subroutine fluxes_centdif

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc = box%n_cell
    inv_dr = 1/box%dr

    ! Add the fluxes to the cells
    box%cc(1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, i_phi) + dt(1) * ( &
         (box%fx(1:nc, :, i_phi) - box%fx(2:nc+1, :, i_phi)) * inv_dr + &
         (box%fy(:, 1:nc, i_phi) - box%fy(:, 2:nc+1, i_phi)) * inv_dr)
  end subroutine update_solution

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    real(dp)                 :: diff
    integer                  :: nc, lvl

    nc   = boxes(id)%n_cell
    lvl  = boxes(id)%lvl

    ! Use a curvature monitor for the refinement criterion
    diff = maxval(abs( &
         boxes(id)%cc(2:nc+1, 1:nc, i_phi) + &
         boxes(id)%cc(0:nc-1, 1:nc, i_phi) + &
         boxes(id)%cc(1:nc, 2:nc+1, i_phi) + &
         boxes(id)%cc(1:nc, 0:nc-1, i_phi) - 4 * &
         boxes(id)%cc(1:nc, 1:nc, i_phi)))

    if (diff > 5.0e-3_dp .and. lvl < 6) then
       ref_flags(id) = a5_do_ref
    else if (diff < 0.05_dp * 5.0e-3_dp) then
       ref_flags(id) = a5_rm_ref
    end if
  end subroutine set_ref_flags

  subroutine prolong_to_new_children(tree, ref_info)
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_prolong1_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box_sides(tree%boxes, id, i_phi, &
               a2_sides_interp, a2_bc_dirichlet)
       end do
    end do
  end subroutine prolong_to_new_children

end program
