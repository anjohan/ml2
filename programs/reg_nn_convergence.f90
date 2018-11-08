program convergence
    use mod_neural_network
    use iso_fortran_env, only: dp => real64
    implicit none

    character(len=*), intent(in), optional :: filename
    integer :: L, num_states, i, j, k, u, t
    real(dp), allocatable :: spins(:,:), couplings(:,:), &
                             energies(:,:), pred_energies(:,:), &
                             test_couplings(:,:), test_energies(:,:)
    class(neural_network), allocatable :: nn
    integer :: max_epochs = 1000
    real(dp) :: tolerance = 0.99d0
    real(dp) :: best_learning_rate, r2
    integer :: fastest_convergence, best_batch_size

    real(dp), allocatable :: lambdas(:), learning_rates(:)
    integer, allocatable :: batch_sizes(:)

    setup: block
        L = 40
        num_states = 1600
        allocate(spins(L, 2*num_states), couplings(L**2, 2*num_states), &
                 energies(1, 2*num_states), &
                 pred_energies(1,num_states))

        call random_number(spins)
        where (spins < 0.5d0)
            spins = -1
        elsewhere
            spins = 1
        end where

        do i = 1, 2*num_states
            associate(s => spins(:,i))
                couplings(:,i) = -spin_coupling_vector(s)
                energies(1,i) = energy(s)
            end associate
        end do

        test_couplings = couplings(:,num_states+1:)
        couplings = couplings(:,:num_states)
        test_energies = energies(:,num_states+1:)
        energies = energies(:,:num_states)

        lambdas = [(10**i, i = -4, 2)]
        num_lambdas = size(lambdas)
        batch_sizes = [1,10,40,100,200,500,1000]
        num_batch_sizes = size(batch_sizes)
        learning_rates = [(10**i, 5*10**i, i = -4, 1)]
        num_learning_rates = size(learning_rates)
    end block setup

    nn = neural_network(L**2, [1], relu())

    open(newunit=u, file="data/reg_nn_convergence.dat", status="replace")

    do i = 1, num_lambdas
        nn%lambda = lambdas(i)
        fastest_convergence = max_epochs
        best_batch_size = batch_sizes(1)
        best_learning_rate = learning_rates(1)
        do j = 1, num_batch_sizes
            do k = 1, num_learning_rates
                call nn%reset_weights()
                do t = 1, ,max_epochs
                    call nn%train(couplings, energies, learning_rates(k), &
                                  1, batch_sizes(j))
                    do j = 1, num_states
                        call nn%predict(test_couplings(:,j), pred_energies(:,j))
                    end do
                    r2 = 1 - sum((pred_energies(1,:) - test_energies(1,:))**2) &
                             / sum((test_energies(1,:) - sum(test_energies(1,:))/num_states)**2)
                    if (r2 >= tolerance .and. t < fastest_convergence) then
                        best_batch_size = batch_sizes(j)
                        best_learning_rate = learning_rates(k)
                end do
            end do
        end do
end program
