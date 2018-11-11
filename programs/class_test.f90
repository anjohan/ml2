program class_test
    use mod_ising, only: read_2d_states
    use mod_activation_functions
    use mod_neural_network
    use iso_fortran_env, only: dp => real64
    implicit none

    class(neural_network), allocatable :: nn1, nn2, nn3

    real(dp), allocatable :: X_train(:,:), X_test(:,:), X_crit(:,:)
    real(dp), allocatable :: y_train(:,:), y_test(:,:), y_crit(:,:)
    integer, allocatable :: y_train_int(:), y_test_int(:), &
                            y_crit_int(:)

    real(dp), allocatable :: train_acc(:,:), test_acc(:,:), &
                             crit_acc(:,:)

    real(dp) :: test_fraction = 0.2d0
    real(dp) :: learning_rates(2) = [0.001d0, 0.01d0], &
                pred1(1), pred2(1), pred3(1)

    integer :: num_spins, num_epochs, i, j, k, u

    call read_2d_states("data/states_10.bin", "data/labels_10.bin", &
                        X_train, X_test, X_crit, &
                        y_train_int, y_test_int, y_crit_int, &
                        test_fraction)
    allocate(y_train(1, size(y_train_int)))
    allocate(y_test(1, size(y_test_int)))
    allocate(y_crit(1, size(y_crit_int)))

    y_train(1,:) = y_train_int(:)
    y_test(1,:) = y_test_int(:)
    y_crit(1,:) = y_crit_int(:)
    num_spins = size(X_train, 1)
    num_epochs = 20

    allocate(train_acc(num_epochs, 3*size(learning_rates)))
    train_acc(:,:) = 0
    allocate(test_acc, crit_acc, source=train_acc)

    nn1 = neural_network(num_spins, [(10, i = 1, 1), 1], relu(), &
                         output_activation=sigmoid(), lambda=0.001d0)
    nn2 = neural_network(num_spins, [10,10,1], relu(), &
                         output_activation=sigmoid(), lambda=0.001d0)
    nn3 = neural_network(num_spins, [100,1], relu(), &
                         output_activation=sigmoid(), lambda=0.001d0)
    do i = 1, size(learning_rates)
        call nn1%reset_weights()
        call nn2%reset_weights()
        call nn3%reset_weights()

        do j = 1, num_epochs
            write(*,*) this_image(), i, j
            call nn1%train(X_train, y_train, &
                           learning_rates(i), 1, 20)
            call nn2%train(X_train, y_train, &
                           learning_rates(i), 1, 20)
            call nn3%train(X_train, y_train, &
                           learning_rates(i), 1, 20)
            do k = 1, size(y_test, 2)
                call nn1%predict(X_test(:,k), pred1)
                call nn2%predict(X_test(:,k), pred2)
                call nn3%predict(X_test(:,k), pred3)
                if ((pred1(1)>0.5d0 .and. y_test_int(k)==1) .or. &
                    (pred1(1)<0.5d0 .and. y_test_int(k)==0)) then
                    test_acc(j, 3*i-2) = test_acc(j, 3*i-2) + 1
                end if
                if ((pred2(1)>0.5d0 .and. y_test_int(k)==1) .or. &
                    (pred2(1)<0.5d0 .and. y_test_int(k)==0)) then
                    test_acc(j, 3*i-1) = test_acc(j, 3*i-1) + 1
                end if
                if ((pred3(1)>0.5d0 .and. y_test_int(k)==1) .or. &
                    (pred3(1)<0.5d0 .and. y_test_int(k)==0)) then
                    test_acc(j, 3*i) = test_acc(j, 3*i) + 1
                end if
            end do
            if (this_image() == 1) write(*,*) test_acc(j,:)/size(y_test_int)
        end do
    end do

    test_acc(:,:) = test_acc(:,:)/size(X_test, 2)
    if (this_image() == 1) then
        open(newunit=u, file="data/class_test.dat", &
             status="replace")
        do i = 1, num_epochs
            write(u, *) i, test_acc(i, :)
        end do
        close(u)
    end if

end program
