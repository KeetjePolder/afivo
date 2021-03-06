!> \example random_refinement_3d.f90

! This example shows how to create an AMR tree, perform random refinement, write
! output files, and how to fill ghost cells.
program random_refinement_3d
  use m_a3_t
  use m_a3_core
  use m_a3_utils
  use m_a3_io
  use m_a3_gc

  implicit none

  type(a3_t)           :: tree
  integer              :: i
  integer, parameter   :: n_boxes_base = 1
  integer              :: ix_list(3, n_boxes_base)
  integer              :: nb_list(6, n_boxes_base)
  integer, parameter   :: n_cell       = 8
  integer, parameter   :: i_phi        = 1
  type(ref_info_t)     :: ref_info
  real(dp)             :: dr
  character(len=40)    :: fname

  ! The cell spacing at the coarsest grid level
  dr = 2 * acos(-1.0_dp) / n_cell ! 2 * pi / n_cell

  ! Initialize tree
  call a3_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       1, &            ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       cc_names = ["phi"])      ! Optional: names of cell-centered variables

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1,1] touches the origin (x,y,z) = (0,0,0)
  ix_list(:, 1) = [1,1,1] ! One box at index 1,1

  ! Set neighbors for box one, here nb means neighbor (direction) and l/h stands
  ! for low/high
  nb_list(a3_neighb_lowx, 1) = 1      ! lower-x neighbor of box 1 is box 1
  nb_list(a3_neighb_highx, 1) = 1      ! higher-x neighbor of box 1 is box 1
  nb_list(a3_neighb_lowy, 1) = 1      ! idem for y-direction
  nb_list(a3_neighb_highy, 1) = 1
  nb_list(a3_neighb_lowz, 1) = 1      ! idem for z-direction
  nb_list(a3_neighb_highz, 1) = 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a3_set_base(tree, ix_list, nb_list)

  ! Set variables on base by using the helper functions a3_loop_box(tree, sub)
  ! and a3_loop_boxes(tree, sub). These functions call the subroutine sub for
  ! each box in the tree, with a slightly different syntax.
  call a3_loop_box(tree, set_init_cond)

  ! Fill ghost cells for phi. The third argument is a subroutine that fills
  ! ghost cells near refinement boundaries, and the fourth argument fill ghost
  ! cells near physical boundaries.
  call a3_gc_tree(tree, i_phi, a3_gc_interp, a3_bc_dirichlet_zero)

  do i = 1, 25
     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement). Variables are the names
     ! given as the third argument.
     write(fname, "(A,I0)") "random_refinement_3d_", i
     call a3_write_vtk(tree, trim(fname), dir="output")

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a3_adjust_refinement(tree, ref_routine, ref_info)
     print *, "# new     boxes", ref_info%n_add
     print *, "# removed boxes", ref_info%n_rm

     ! Newly added boxes should be initialized, which is done in this routine.
     call prolong_to_new_children(tree, ref_info)
  end do

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a3_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:) ! A list of all boxes in the tree
    integer, intent(in)      :: id       ! The index of the current box
    integer, intent(inout)   :: ref_flag
    real(dp)                 :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.5_dp**0.125_dp .and. boxes(id)%lvl < 4) then
       ref_flag = a5_do_ref ! Add refinement
    else
       ref_flag = a5_rm_ref ! Ask to remove this box, which will not always
                            ! happen (see documentation)
    end if
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

             ! Set the values at each cell according to some function
             box%cc(i, j, k, i_phi) = sin(0.5_dp * xyz(1)) * cos(xyz(2))
          end do
       end do
    end do
  end subroutine set_init_cond

  ! Set values on newly added boxes
  subroutine prolong_to_new_children(tree, ref_info)
    use m_a3_prolong
    type(a3_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       ! For each newly added box ...
       do i = 1, size(ref_info%lvls(lvl)%add)
          ! Use prolongation to set the values of variable i_phi
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent

          call a3_prolong2(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          ! After values have been set on this level, fill ghost cells
          call a3_gc_box(tree%boxes, id, i_phi, &
               a3_gc_interp, a3_bc_dirichlet_zero)
       end do
    end do
  end subroutine prolong_to_new_children

end program random_refinement_3d
