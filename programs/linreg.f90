program linreg
    use iso_fortran_env, only: dp => real64
    use mod_ols
    use mod_ridge
    use mod_ising_basis
    use mod_ising
    use mod_bootstrap
    implicit none

    class(bootstrapper), allocatable :: bs
    class(regressor), allocatable :: fitter
    class(coupling), allocatable :: basis(:)

    integer :: L, num_states, num_lambdas, num_bootstraps, i, u, k, u_bivar
    real(dp), allocatable :: spins(:,:), couplings(:,:), energies(:), J(:,:), &
                             lambdas(:), test_couplings(:,:), test_energies(:), tmp(:)
    real(dp) :: mse, r2
    character(len=1024) :: str

    num_bootstraps = 20
    L = 40

    read(*,*) num_states
    write(str,"(i0)") num_states

    call create_ising_basis(basis, L)

    allocate(spins(L, 2*num_states), couplings(2*num_states, L**2), &
             energies(2*num_states), tmp(num_states))

    call random_number(spins)
    where (spins < 0.5d0)
        spins = -1
    elsewhere
        spins = 1
    end where

    do i = 1, 2*num_states
        associate(s => spins(:, i))
            couplings(i, :) = spin_coupling_vector(s)
            energies(i) = energy(s)
        end associate
    end do

    test_couplings = couplings(num_states+1:, :)
    test_energies = energies(num_states+1:)
    couplings = couplings(1:num_states, :)
    energies = energies(1:num_states)

    fitter = ols(basis)
    call fitter%fit(couplings, energies)
    J = reshape(fitter%beta, [L, L])

    open(newunit=u, file="data/J_ols_" // trim(str) // ".dat", status="replace")
    write(u, "(i0,/)") L
    do i = 1, L
        write(u, "(*(f0.6,:,x))") J(:, i)
    end do
    close(u)

    ! ==== Ridge ==== !
    lambdas = 10.0d0**[(i, i = -4, 4)]
    num_lambdas = size(lambdas)

    open(newunit=u, file="data/J_Ridge_" // trim(str) // ".dat", status="replace")
    write(u, *) L

    open(newunit=u_bivar, file="data/mse_Ridge_" // trim(str) // ".dat", status="replace")

    do i = 1, num_lambdas
        fitter = ridge(lambdas(i), basis)
        call fitter%fit(couplings, energies)
        J = reshape(fitter%beta, [L, L])


        write(u, *)
        do k = 1, L
            write(u, "(*(f0.6,:,x))") J(:, k)
        end do

        bs = bootstrapper(fitter)
        call bs%bootstrap(couplings, energies, num_bootstraps, 0.2d0)

        fitter%beta = bs%mean_beta
        call fitter%predict(test_couplings, tmp, test_energies, mse, r2)

        write(u_bivar, *) lambdas(i), r2, bs%mean_MSE, bs%bias, bs%variance
    end do

    close(u)
    close(u_bivar)

    open(newunit=u, file="data/Ridge_lambdas.dat", status="replace")
    write(u, "(*(es8.2,:,/))") lambdas
    close(u)
end program
