program logreg
    use mod_ising
    use mod_utilities, only: shuffle
    use mod_binary_logreg
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), allocatable :: X_train(:,:), X_test(:,:), X_crit(:,:)
    integer, allocatable :: y_test(:), y_train(:), y_pred(:), y_crit(:)
    integer :: num_spins, num_crit, num_noncrit
    integer :: u, p
    real(dp) :: test_fraction
    integer :: N_test

    test_fraction = 0.2d0

    write(*,*) "Using", num_images(), " images"

    call read_2d_states("data/states.bin", "data/labels.bin", &
                        X_train, X_test, X_crit, y_train, y_test, y_crit, &
                        test_fraction, add_intercept=.true., T=.true.)


    simulations: block
        class(binary_logreg), allocatable :: fitter

        real(dp), allocatable :: lambdas(:), momentums(:), learning_rates(:), &
                                 training_accuracies(:,:,:), &
                                 test_accuracies(:,:,:), &
                                 crit_accuracies(:,:,:)
        integer :: num_lambdas, num_momentums, num_learning_rates, i, j, k, iteration
        integer, allocatable :: training_pred(:), test_pred(:), crit_pred(:)

        lambdas = [0.0d0, (10.0d0**i, i = -5,2)]
        momentums = [(0.2d0*i, i = 0, 5)]
        learning_rates = [(10.0d0**i, i = -4, 0)]

        num_lambdas = size(lambdas)
        num_momentums = size(momentums)
        num_learning_rates = size(learning_rates)

        allocate(training_accuracies(num_lambdas, num_learning_rates, num_momentums))
        training_accuracies(:,:,:) = 0
        allocate(test_accuracies, crit_accuracies, source=training_accuracies)

        allocate(training_pred, mold=y_train)
        allocate(test_pred, mold=y_test)
        allocate(crit_pred, mold=y_crit)

        fitter = binary_logreg(learning_rate=1.0d0, batch_size=32)

        iteration = -1

        !!$omp parallel do collapse(3) default(shared) &
        !!$omp& private(i,j,k,test_pred,training_pred,crit_pred) firstprivate(fitter)
        do k = 1, num_momentums
            do j = 1, num_learning_rates
                do i = 1, num_lambdas
                    iteration = iteration + 1
                    if (mod(iteration, num_images()) /= this_image()-1) cycle

                    write(*, "('Image: ', i0, 3(', ', i0, ' of ', i0))") this_image(), &
                        k, num_momentums, j, num_learning_rates, i, num_lambdas

                    fitter%lambda = lambdas(i)
                    fitter%learning_rate = learning_rates(j)
                    fitter%momentum = momentums(k)

                    call fitter%fit(X_train, y_train)
                    call fitter%predict(X_train, training_pred)
                    call fitter%predict(X_test, test_pred)
                    call fitter%predict(X_crit, crit_pred)

                    training_accuracies(i,j,k) = count(y_train == training_pred) &
                                                 / (1.0d0*size(y_train))
                    test_accuracies(i,j,k) = count(y_test == test_pred) &
                                                 / (1.0d0*size(y_test))
                    crit_accuracies(i,j,k) = count(y_crit == crit_pred) &
                                                 / (1.0d0*size(y_crit))
                    write(*,*) training_accuracies(i,j,k), test_accuracies(i,j,k), &
                                crit_accuracies(i,j,k)
                end do
            end do
        end do
        !!$omp end parallel do
        call co_sum(training_accuracies)
        call co_sum(test_accuracies)
        call co_sum(crit_accuracies)

        if (this_image() == 1) then
        post_analysis: block
            integer :: best_indices(2)
            character(len=1024) :: filename
            real(dp), allocatable :: tot_accuracy(:,:,:)

            tot_accuracy = test_accuracies + crit_accuracies

            open(newunit=u, file="data/logreg_table.dat", status="replace")
            write(u, "(a)") "lambda {Training accuracy} {Test accuracy}" &
                            // " {Critical accuracy} {Learning rate} {Momentum}"

            do i = 1, num_lambdas
                best_indices = maxloc(tot_accuracy(i,:,:))
                j = best_indices(1); k = best_indices(2)

                write(u, "(es10.0,x,f0.2,x,f0.2,x,f0.2,x,es10.0,x,f0.2)") &
                    lambdas(i), training_accuracies(i,j,k), test_accuracies(i,j,k), &
                    crit_accuracies(i,j,k), learning_rates(j), momentums(k)
            end do
            close(u)
        end block post_analysis
        end if

    end block simulations

end program
