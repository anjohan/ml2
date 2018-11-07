program reg_nn
    use mod_ising
    use mod_neural_network
    use iso_fortran_env, only: dp => real64

    integer :: L, num_states, i, j
    real(dp), allocatable :: spins(:,:), couplings(:,:), energies(:,:), pred_energies(:,:)
    class(neural_network), allocatable :: nn
    real(dp) :: cost, r2


    L = 40
    num_states = 1000
    allocate(spins(L, num_states), couplings(L**2, num_states), energies(1, num_states), &
             pred_energies(1,num_states))

    call random_number(spins)
    where (spins < 0.5d0)
        spins = -1
    elsewhere
        spins = 1
    end where

    do i = 1, num_states
        associate(s => spins(:,i))
            couplings(:,i) = -spin_coupling_vector(s)
            energies(1,i) = energy(s)
        end associate
    end do

    nn = neural_network(L, [100,100,1], relu())

    do i = 1, 100
        call nn%train(spins, energies, 0.1d0, 1, num_states/100)
        cost = nn%cost_func(spins, energies, l2norm)
        do j = 1, num_states
            call nn%predict(spins(:,j), pred_energies(:,j))
        end do
        r2 = 1 - sum((pred_energies(1,:) - energies(1,:))**2) &
             / sum((energies(1,:) - sum(energies(1,:))/N)**2)
        write(*,*) cost, r2
    end do
end program
