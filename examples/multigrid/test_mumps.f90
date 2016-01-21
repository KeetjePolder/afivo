program test_mumps
  use m_mumps

  implicit none

  call simple_test()
  call tree_test()

contains

  subroutine simple_test()
    type(DMUMPS_STRUC) :: mp

    call mumps_init(mp, 5, 12)

    mp%RHS = [20.0_dp, 24.0_dp, 9.0_dp, 6.0_dp, 13.0_dp]
    mp%IRN = [1, 2, 4, 5, 2, 1, 5, 3, 2, 3, 1, 3]
    mp%JCN = [2, 3, 3, 5, 1, 1, 2, 4, 5, 2, 3, 3]
    mp%A = [3.0_dp, -3.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, &
         4.0_dp, 2.0_dp, 6.0_dp, -1.0_dp, 4.0_dp, 1.0_dp]

    call mumps_factor(mp)
    call mumps_solve(mp)
    print *, "Solution is ", mp%RHS
    print *, "Should be 1,2,3,4,5"
    call mumps_free(mp)
  end subroutine simple_test

  subroutine tree_test()
    use m_a2_t
    use m_a2_core
    use m_a2_mg
    use m_a2_gc

    integer, parameter :: box_size = 2
    integer, parameter :: i_phi = 1, i_tmp = 2
    integer, parameter :: i_rhs = 3
    real(dp), parameter :: dr_base = 1.0
    integer            :: ix_list(2, 1)
    integer            :: nb_list(4, 1)

    type(a2_t)         :: tree
    type(mg2_t)        :: mg
    type(DMUMPS_STRUC) :: mp

    call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, dr=dr_base)

    ix_list(:, 1) = [1,1]         ! Set index of boxnn
    nb_list(:, 1) = -1            ! Dirichlet zero -> -1
    call a2_set_base(tree, ix_list, nb_list)

    call a2_print_info(tree)

    mg%i_phi        = i_phi
    mg%i_tmp        = i_tmp
    mg%i_rhs        = i_rhs
    mg%sides_bc     => a2_bc_dirichlet_zero

    call mg2_init_mg(mg)
    call mumps_factor_2d(mp, tree, mg)
  end subroutine tree_test
end program test_mumps
