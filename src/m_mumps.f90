module m_mumps

  implicit none
  public

  include 'mpif.h'
  include 'dmumps_struc.h'

  integer, parameter :: dp = kind(0.0d0)

contains

  subroutine mumps_init(mp, n_vars, n_nonzeros)
    type(DMUMPS_STRUC), intent(inout) :: mp
    integer, intent(in) :: n_vars, n_nonzeros

    mp%COMM = MPI_COMM_WORLD
    mp%JOB = -1
    mp%SYM = 0
    mp%PAR = 1
    call DMUMPS(mp)

    if (mp%INFOG(1) < 0) then
       print *, "MUMPS initialize error: ", mp%INFOG(1:2)
       stop
    end if

    mp%N = n_vars
    mp%NZ = n_nonzeros

    allocate(mp%IRN(n_nonzeros))
    allocate(mp%JCN(n_nonzeros))
    allocate(mp%A(n_nonzeros))
    allocate(mp%RHS(n_vars))

    mp%ICNTL(4) = 1             ! Only enable errors
  end subroutine mumps_init

  subroutine mumps_factor(mp)
    type(DMUMPS_STRUC), intent(inout) :: mp

    mp%JOB = 4
    call DMUMPS(mp)

    if (mp%INFOG(1) < 0) then
       print *, "MUMPS factorization error: ", mp%INFOG(1:2)
       stop
    end if
  end subroutine mumps_factor

  subroutine mumps_solve(mp)
    type(DMUMPS_STRUC), intent(inout) :: mp

    mp%JOB = 3
    call DMUMPS(mp)

    if (mp%INFOG(1) < 0) then
       print *, "MUMPS solve error: ", mp%INFOG(1:2)
       stop
    end if
  end subroutine mumps_solve

  subroutine mumps_free(mp)
    type(DMUMPS_STRUC), intent(inout) :: mp

    mp%JOB = -2
    call DMUMPS(mp)

    if (mp%INFOG(1) < 0) then
       print *, "MUMPS solve error: ", mp%INFOG(1:2)
       stop
    end if

    deallocate(mp%IRN)
    deallocate(mp%JCN)
    deallocate(mp%A)
    deallocate(mp%RHS)
  end subroutine mumps_free

  subroutine mumps_factor_2d(mp, tree, mg)
    use m_a2_t
    use m_a2_mg

    type(DMUMPS_STRUC), intent(inout) :: mp
    type(a2_t), intent(in)            :: tree
    type(mg2_t), intent(in)           :: mg

    if (.not. associated(mg%box_op, mg2_box_lpl)) stop "only implemented for laplacian"

    nc = tree%n_cell

    ! Count number of boxes at base level
    n_boxes = size(tree%lvls(1)%ids)

    ! Number of unknowns
    n_vars = n_boxes * nc**2

    ! Number of non-zeros in the matrix

    ! Store only half (symmetric matrix)
  end subroutine mumps_factor_2d

end module m_mumps
