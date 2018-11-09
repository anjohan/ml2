program convergence
    use mod_ising
    use mod_neural_network
    use iso_fortran_env, only: dp => real64
    implicit none

    integer :: L, num_states, i, j, k, u, t, m
    real(dp), allocatable :: spins(:,:)[:], couplings(:,:), &
                             energies(:,:), pred_energies(:,:), &
                             test_couplings(:,:), test_energies(:,:)
    class(neural_network), allocatable :: nn
    integer :: max_epochs = 100
    real(dp) :: tolerance = 0.95d0
    real(dp) :: best_learning_rate, r2
    integer :: fastest_convergence, best_batch_size
    integer :: num_batch_sizes, num_lambdas, num_learning_rates

    real(dp), allocatable :: lambdas(:), learning_rates(:)
    integer, allocatable :: batch_sizes(:)

    setup: block
        L = 40
        num_states = 1600
        allocate(spins(L, 2*num_states)[*], couplings(L**2, 2*num_states), &
                 energies(1, 2*num_states), &
                 pred_energies(1,num_states))

        call random_number(spins)
        sync all
        spins(:,:) = spins(:,:)[1]
        where (spins < 0.5d0)
            spins = -1
        elsewhere
            spins = 1
        end where

        do i = 1, 2*num_states
            couplings(:,i) = -spin_coupling_vector(spins(:,i))
            energies(1,i) = energy(spins(:,i))
        end do

        test_couplings = couplings(:,num_states+1:)
        couplings = couplings(:,:num_states)
        test_energies = energies(:,num_states+1:)
        energies = energies(:,:num_states)

        lambdas = [0.0d0, (10.0d0**i, i = -5, -1,3)]
        num_lambdas = size(lambdas)
        batch_sizes = [1,10,40,100] !,200,500,1000]
        num_batch_sizes = size(batch_sizes)
        learning_rates = [(10.0d0**i, 5*10.0d0**i, i = -5, -1, 2)]
        num_learning_rates = size(learning_rates)
    end block setup

    nn = neural_network(L**2, [1], relu())

    if (this_image() == 1) then
        open(newunit=u, file="data/reg_nn_convergence.dat", status="replace")
        write(u,*) "lambda {Learning rate} {Batch size} {Epochs for convergence}"
    end if

    do i = 1, num_lambdas
        nn%lambda = lambdas(i)
        fastest_convergence = max_epochs
        best_batch_size = batch_sizes(1)
        best_learning_rate = learning_rates(1)
        do j = 1, num_batch_sizes
            params: do k = 1, num_learning_rates
                !if (this_image() == 1) write(*,*) fastest_convergence
                call nn%reset_weights()
                do t = 1,max_epochs
                    call nn%train(couplings, energies, learning_rates(k), &
                                  1, batch_sizes(j))
                    do m = 1, num_states
                        call nn%predict(test_couplings(:,m), pred_energies(:,m))
                    end do
                    r2 = 1 - sum((pred_energies(1,:) - test_energies(1,:))**2) &
                             / sum((test_energies(1,:) - sum(test_energies(1,:))/num_states)**2)
                    !write(*,*) this_image(),i,j,k,t, r2
                    call co_max(r2)
                    if (r2 >= tolerance) then
                        if (t < fastest_convergence) then
                            best_batch_size = batch_sizes(j)
                            best_learning_rate = learning_rates(k)
                            fastest_convergence = t
                            if (this_image() == 1) write(*,*) t, i, batch_sizes(j), learning_rates(k), r2
                        end if
                        cycle params
                    end if
                end do
            end do params
        end do
        if (this_image() == 1) write(u, "(es10.1,es10.1,x,i0,x,i0)") lambdas(i), &
                            best_learning_rate, best_batch_size, fastest_convergence
    end do
    if (this_image() == 1) close(u)
end program
