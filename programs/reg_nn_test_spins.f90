program reg_nn_spins
    use mod_ising
    use mod_neural_network
    use iso_fortran_env, only: dp => real64
    implicit none

    integer :: L, num_states, i, j, k, num_epochs, u
    class(neural_network), allocatable :: nn

    real(dp), allocatable :: learning_rates(:)
    real(dp), allocatable :: spins(:,:)[:], energies(:,:), pred(:,:)
    real(dp), allocatable :: test_spins(:,:)[:], test_energies(:,:), r2s(:,:)

    L = 40
    num_states = 10000
    num_epochs = 100
    allocate(spins(L, num_states)[*], energies(1, num_states), pred(1, num_states))
    allocate(test_spins(L, num_states)[*], test_energies(1, num_states))

    call random_number(spins); call random_number(test_spins)

    where (spins > 0.5)
        spins = 1
    elsewhere
        spins = -1
    end where
    where (test_spins > 0.5)
        test_spins = 1
    elsewhere
        test_spins = -1
    end where


    sync all
    spins(:,:) = spins(:,:)[1]; test_spins(:,:) = test_spins(:,:)[1]
    do i = 1, num_states
        energies(1, i) = energy(spins(:, i))
        test_energies(1, i) = energy(test_spins(:, i))
    end do

    learning_rates = [0.001d0, 0.01d0, 0.04d0]
    allocate(r2s(num_epochs, size(learning_rates)))

    nn = neural_network(L, [50,50,50,1], relu(), lambda=0.001d0)
    do i = 1, size(learning_rates)
        call nn%reset_weights()
        do j = 1, num_epochs
            write(*,*) i, j
            call nn%train(spins, energies, learning_rates(i), 1, 32)
            do k = 1, num_states
                call nn%predict(test_spins(:, k), pred(:, k))
            end do
            r2s(j, i) = 1 - sum((pred(1,:) - test_energies(1,:))**2) &
                      / sum((test_energies(1,:) - sum(test_energies(1,:))/num_states)**2)
            if (this_image() == 1) write(*,*) j, r2s(j,i)
        end do
    end do

    if (this_image() == 1) then
        open(newunit=u, file="data/reg_nn_test_spins.dat", status="replace")
        do i = 1, num_epochs
            write(u, *) i, r2s(i, :)
        end do
        close(u)
    end if

end program
