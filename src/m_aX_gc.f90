! This module contains routines related to the filling of ghost cells. Note that
! corner ghost cells are not used in Afivo.
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_gc
  use m_a$D_t

  implicit none
  private

  public :: a$D_gc_tree
  public :: a$D_gc_box
  public :: a$D_bc_dirichlet_zero
  public :: a$D_bc_neumann_zero
  public :: a$D_gc_interp
  public :: a$D_gc_prolong0
  public :: a$D_gc_interp_lim
  public :: a$D_gc2_box
  public :: a$D_gc2_prolong1
  public :: a$D_bc2_neumann_zero

contains

  !> Fill ghost cells for variables iv on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine a$D_gc_tree(tree, iv, subr_rb, subr_bc)
    type(a$D_t), intent(inout) :: tree !< Tree to fill ghost cells on
    integer, intent(in)       :: iv !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)     :: subr_rb !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)     :: subr_bc    !< Procedure called at physical boundaries
    integer                   :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a$D_gc_box(tree%boxes, id, iv, subr_rb, subr_bc)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_gc_tree

  !> Fill ghost cells for variables iv on the sides of a box, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine a$D_gc_box(boxes, id, iv, subr_rb, subr_bc)
    type(box$D_t), intent(inout)  :: boxes(:)              !< List of all the boxes
    integer, intent(in)          :: id                    !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv                    !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)       :: subr_rb               !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)       :: subr_bc               !< Procedure called at physical boundaries
    integer                      :: nb, nb_id, bc_type

    do nb = 1, a$D_num_neighbors
       nb_id = boxes(id)%neighbors(nb)
       if (nb_id > a5_no_box) then
          call sides_from_nb(boxes(id), boxes(nb_id), nb, iv)
       else if (nb_id == a5_no_box) then
          call subr_rb(boxes, id, nb, iv)
       else
          call subr_bc(boxes(id), nb, iv, bc_type)
          call bc_to_gc(boxes(id), nb, iv, bc_type)
       end if
    end do
  end subroutine a$D_gc_box

  subroutine bc_to_gc(box, nb, iv, bc_type)
    type(box$D_t), intent(inout)  :: box
    integer, intent(in)          :: iv                    !< Variable to fill
    integer, intent(in)          :: nb                    !< Neighbor direction
    integer, intent(in)          :: bc_type               !< Type of b.c.
    real(dp)                     :: c1, c2
    integer                      :: nc

    nc = box%n_cell

    ! If we call the interior point x, and the ghost point y, then a Dirichlet
    ! boundary value b can be imposed as:
    ! y = -x + 2*b
    ! A Neumann b.c. can be imposed as y = x +/- dx * b
    ! Below, we set coefficients c1 and c2 to handle both cases
    select case (bc_type)
    case (a5_bc_dirichlet)
       c1 = -1
       c2 = 2
    case (a5_bc_neumann)
       c1 = 1
       c2 = box%dr * a$D_neighb_high_pm(nb) ! This gives a + or - sign
    case default
       stop "fill_bc: unknown boundary condition"
    end select

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       box%cc(0, 1:nc, iv) = &
            c2 * box%cc(0, 1:nc, iv) + c1 * box%cc(1, 1:nc, iv)
    case (a2_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = &
            c2 * box%cc(nc+1, 1:nc, iv) + c1 * box%cc(nc, 1:nc, iv)
    case (a2_neighb_lowy)
       box%cc(1:nc, 0, iv) = &
            c2 * box%cc(1:nc, 0, iv) + c1 * box%cc(1:nc, 1, iv)
    case (a2_neighb_highy)
       box%cc(1:nc, nc+1, iv) = &
            c2 * box%cc(1:nc, nc+1, iv) + c1 * box%cc(1:nc, nc, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv) = &
            c2 * box%cc(0, 1:nc, 1:nc, iv) + c1 * box%cc(1, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = &
            c2 * box%cc(nc+1, 1:nc, 1:nc, iv) + c1 * box%cc(nc, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv) = &
            c2 * box%cc(1:nc, 0, 1:nc, iv) + c1 * box%cc(1:nc, 1, 1:nc, iv)
    case (a3_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = &
            c2 * box%cc(1:nc, nc+1, 1:nc, iv) + c1 * box%cc(1:nc, nc, 1:nc, iv)
    case (a3_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv) = &
            c2 * box%cc(1:nc, 1:nc, 0, iv) + c1 * box%cc(1:nc, 1:nc, 1, iv)
    case (a3_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = &
            c2 * box%cc(1:nc, 1:nc, nc+1, iv) + c1 * box%cc(1:nc, 1:nc, nc, iv)
#endif
    end select
  end subroutine bc_to_gc

  !> Partial prolongation to the ghost cells of box id from parent
  subroutine a$D_gc_prolong0(boxes, id, nb, iv)
    use m_a$D_prolong, only: a$D_prolong0
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in)           :: nb       !< Neighbor to get data from
    integer                       :: p_id, nb_dim, lo($D), hi($D)

    p_id       = boxes(id)%parent
    nb_dim     = a$D_neighb_dim(nb)
    lo(:)      = 1
    hi(:)      = boxes(id)%n_cell
    lo(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell+1)
    hi(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell+1)

    call a$D_prolong0(boxes(p_id), boxes(id), iv, low=lo, high=hi)
  end subroutine a$D_gc_prolong0

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine a$D_gc_interp(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_neighb_dim(nb))
#if $D == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          boxes(id)%cc(ix, j, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(ix_c, j_c1, iv) + &
               sixth * boxes(p_nb_id)%cc(ix_c, j_c2, iv) + &
               third * boxes(id)%cc(ix_f, j, iv)
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, ix, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, ix_c, iv) + &
               sixth * boxes(p_nb_id)%cc(i_c2, ix_c, iv) + &
               third * boxes(id)%cc(i, ix_f, iv)
       end do
#elif $D==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             boxes(id)%cc(ix, j, k, iv) = &
                  third * boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv) + &
                  third * boxes(id)%cc(ix_f, j, k, iv)
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, ix, k, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv) + &
                  third * boxes(id)%cc(i, ix_f, k, iv)
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, ix, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv) + &
                  third * boxes(id)%cc(i, j, ix_f, iv)
          end do
       end do
#endif
    end select

  end subroutine a$D_gc_interp

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries. The ghost values are less than twice the coarse
  !> values.
  subroutine a$D_gc_interp_lim(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp)                    :: c1, c2
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
    real(dp)                    :: c3
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_neighb_dim(nb))
#if $D == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(ix_c, j_c1, iv)
          c2 = boxes(p_nb_id)%cc(ix_c, j_c2, iv)
          boxes(id)%cc(ix, j, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(ix_f, j, iv)
          if (boxes(id)%cc(ix, j, iv) > 2 * c1) boxes(id)%cc(ix, j, iv) = 2 * c1
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(i_c1, ix_c, iv)
          c2 = boxes(p_nb_id)%cc(i_c2, ix_c, iv)
          boxes(id)%cc(i, ix, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(i, ix_f, iv)
          if (boxes(id)%cc(i, ix, iv) > 2 * c1) boxes(id)%cc(i, ix, iv) = 2 * c1
       end do
#elif $D==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv)
             boxes(id)%cc(ix, j, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(ix_f, j, k, iv)
             if (boxes(id)%cc(ix, j, k, iv) > 2 * c1) &
                  boxes(id)%cc(ix, j, k, iv) = 2 * c1
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv)
             boxes(id)%cc(i, ix, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, ix_f, k, iv)
             if (boxes(id)%cc(i, ix, k, iv) > 2 * c1) &
                  boxes(id)%cc(i, ix, k, iv) = 2 * c1
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv)
             c2 = boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv)
             c3 = boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv)
             boxes(id)%cc(i, j, ix, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, j, ix_f, iv)
             if (boxes(id)%cc(i, j, ix, iv) > 2 * c1) &
                  boxes(id)%cc(i, j, ix, iv) = 2 * c1
          end do
       end do
#endif
    end select

  end subroutine a$D_gc_interp_lim

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc_neumann_zero(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb, iv
    integer, intent(out)        :: bc_type

    bc_type = a5_bc_neumann
    call a$D_set_box_gc(box, nb, iv, 0.0_dp)
  end subroutine a$D_bc_neumann_zero

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb, iv
    integer, intent(out)        :: bc_type

    bc_type = a5_bc_dirichlet
    call a$D_set_box_gc(box, nb, iv, 0.0_dp)
  end subroutine a$D_bc_dirichlet_zero

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides_from_nb(box, box_nb, nb, iv)
    type(box$D_t), intent(inout) :: box    !< Box on which to fill ghost cells
    type(box$D_t), intent(in)    :: box_nb !< Neighbouring box
    integer, intent(in)         :: nb        !< Ghost cell / neighbor direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       box%cc(0, 1:nc, iv)    = box_nb%cc(nc, 1:nc, iv)
    case (a2_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = box_nb%cc(1, 1:nc, iv)
    case (a2_neighb_lowy)
       box%cc(1:nc, 0, iv)    = box_nb%cc(1:nc, nc, iv)
    case (a2_neighb_highy)
       box%cc(1:nc, nc+1, iv) = box_nb%cc(1:nc, 1, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv)    = box_nb%cc(nc, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = box_nb%cc(1, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv)    = box_nb%cc(1:nc, nc, 1:nc, iv)
    case (a3_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = box_nb%cc(1:nc, 1, 1:nc, iv)
    case (a3_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv)    = box_nb%cc(1:nc, 1:nc, nc, iv)
    case (a3_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = box_nb%cc(1:nc, 1:nc, 1, iv)
#endif
    end select
  end subroutine sides_from_nb

  !> Get a second layer of ghost cell data (the 'normal' routines give just one
  !> layer of ghost cells). Use subr_rb > on refinement boundaries and subr_bc
  !> on physical boundaries.
  subroutine a$D_gc2_box(boxes, id, iv, subr_rb, subr_bc, gc_data, nc)
    type(box$D_t), intent(inout) :: boxes(:)        !< List of all the boxes
    integer, intent(in)          :: id              !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv              !< Variable for which ghost cells are set
    procedure(a$D_subr_egc)      :: subr_rb         !< Procedure called at refinement boundaries
    procedure(a$D_subr_egc)      :: subr_bc         !< Procedure called at physical boundaries
    integer, intent(in)          :: nc              !< box%n_cell
#if $D   == 2
    real(dp), intent(out)        :: gc_data(nc, 2*$D)     !< The requested ghost cells
#elif $D == 3
    real(dp), intent(out)        :: gc_data(nc, nc, 2*$D) !< The requested ghost cells
#endif
    integer                      :: nb, nb_id

    do nb = 1, a$D_num_neighbors
       nb_id = boxes(id)%neighbors(nb)
       if (nb_id > a5_no_box) then
#if $D == 2
          call sides2_from_nb(boxes(nb_id), nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call sides2_from_nb(boxes(nb_id), nb, iv, gc_data(:, :, nb), nc)
#endif
       else if (nb_id == a5_no_box) then
#if $D == 2
          call subr_rb(boxes, id, nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call subr_rb(boxes, id, nb, iv, gc_data(:, :, nb), nc)
#endif
       else
#if $D == 2
          call subr_bc(boxes, id, nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call subr_bc(boxes, id, nb, iv, gc_data(:, :, nb), nc)
#endif
       end if
    end do
  end subroutine a$D_gc2_box

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides2_from_nb(box_nb, nb, iv, gc_side, nc)
    type(box$D_t), intent(in) :: box_nb !< Neighbouring box
    integer, intent(in)       :: nb     !< Ghost cell / neighbor direction
    integer, intent(in)       :: iv     !< Ghost cell variable
    integer, intent(in)       :: nc
#if $D == 2
    real(dp), intent(out)     :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)     :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, iv)
    case (a2_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, iv)
    case (a2_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, iv)
    case (a2_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, 1:nc, iv)
    case (a3_neighb_lowz)
       gc_side = box_nb%cc(1:nc, 1:nc, nc-1, iv)
    case (a3_neighb_highz)
       gc_side = box_nb%cc(1:nc, 1:nc, 2, iv)
#endif
    end select
  end subroutine sides2_from_nb

  !> Linear interpolation (using data from neighbor) to fill ghost cells
  subroutine a$D_gc2_prolong1(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer, intent(in)         :: nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif
    integer                     :: ix, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)
    ix        = a$D_neighb_high_01(nb) * (nc+3) - 1 ! -1 or nc+2

    select case (a$D_neighb_dim(nb))
#if $D == 2
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          gc_side(j) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          gc_side(i) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
#elif $D==3
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             gc_side(j, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (3)
       k_c1 = ix_offset(3) + ishft(ix+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, j) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
#endif
    end select

  end subroutine a$D_gc2_prolong1

  ! This fills second ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc2_neumann_zero(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, iv)
    case (a2_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, iv)
    case (a2_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, iv)
    case (a2_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (a3_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_neighb_lowz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (a3_neighb_highz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine a$D_bc2_neumann_zero

end module m_a$D_gc
