module mod_find_r2s
    use mod_ising
    use mod_neural_network
    use iso_fortran_env, only: dp => real64
    implicit none

    contains
        function find_r2s(learning_rate, filename) result(r2s)
            real(dp), intent(in) :: learning_rate
            character(len=*), intent(in), optional :: filename
            integer, parameter :: num_epochs = 100
            real(dp) :: r2s(num_epochs)
            integer :: L, num_states, i, j, u
            real(dp), allocatable :: spins(:,:)[:], couplings(:,:), &
                                     energies(:,:), pred_energies(:,:), &
                                     test_couplings(:,:), test_energies(:,:)
            class(neural_network), allocatable :: nn

            L = 40
            num_states = 3000
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

            nn = neural_network(L**2, [1], relu())

            do i = 1, 100
                call nn%train(couplings, energies, learning_rate, 1, num_states/10)
                do j = 1, num_states
                    call nn%predict(test_couplings(:,j), pred_energies(:,j))
                end do
                r2s(i) = 1 - sum((pred_energies(1,:) - test_energies(1,:))**2) &
                         / sum((test_energies(1,:) - sum(test_energies(1,:))/num_states)**2)
            end do

            if (this_image() == 1) then
                if (.not. present(filename)) return
                open(newunit=u, file=filename, status="replace")
                write(u,*) L
                write(u,*)
                do i = 1, L
                    write(u, "(*(f0.6,:,x))") nn%layers(1)%W(:,(i-1)*L+1:i*L)
                end do
                close(u)
            end if
        end function
end module

program reg_nn
    use mod_find_r2s
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp) :: learning_rates(3) = [0.001, 0.01, 0.04]
    real(dp) :: r2s(100,3)
    integer :: i, u

    r2s(:,1) = find_r2s(learning_rates(1), "data/J_nn.dat")
    r2s(:,2) = find_r2s(learning_rates(2))
    r2s(:,3) = find_r2s(learning_rates(3))

    if (this_image() == 1) then
        open(newunit=u, file="data/reg_nn_test_couplings.dat", status="replace")
        do i = 1, 100
            write(u,*) i, r2s(i,:)
        end do
        close(u)
    end if
end program
