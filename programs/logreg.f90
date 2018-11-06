program logreg
    use mod_ising
    use mod_utilities, only: shuffle
    use mod_binary_logreg
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, allocatable :: X_crit_tmp(:,:), X_noncrit(:,:)
    integer, allocatable :: y_crit(:), y_noncrit(:)
    real(dp), allocatable :: X_train(:,:), X_test(:,:), X_crit(:,:)
    integer, allocatable :: y_test(:), y_train(:), y_pred(:)
    integer :: num_spins, num_ordered, num_crit, num_disordered, num_noncrit
    integer :: u, p
    real(dp) :: test_fraction
    integer :: N_test

    test_fraction = 0.2d0

    write(*,*) "Using", num_images(), " images"

    open(newunit=u, file="data/states.bin", action="read", access="stream")
    read(u) num_spins, num_ordered, num_crit, num_disordered
    write(*,*) num_spins, num_ordered, num_crit, num_disordered

    p = num_spins+1
    num_noncrit = num_ordered + num_disordered
    N_test = nint(test_fraction*num_noncrit)

    allocate(X_crit_tmp(num_spins,num_crit), &
             X_noncrit(num_spins,num_noncrit), y_noncrit(num_noncrit))



    allocate(y_crit(num_crit), X_crit(num_crit, p))


        read(u) X_noncrit(:,:num_ordered)
        read(u) X_crit_tmp(:,:)
        read(u) X_noncrit(:,num_ordered+1:)
        close(u)

        open(newunit=u, file="data/labels.bin", access="stream", status="old")
        read(u) y_noncrit(1:num_ordered), y_crit(:), y_noncrit(num_ordered+1:)
        close(u)

        X_noncrit = transpose(X_noncrit)
        X_crit_tmp = transpose(X_crit_tmp)
        X_crit(:,1) = 1
        X_crit(:,2:) = X_crit_tmp(:,:)
        deallocate(X_crit_tmp)
        write(*,*) shape(X_noncrit), shape(y_noncrit)
        write(*,*) shape(X_crit), shape(y_crit)

        write(*,"(*(i0,:,x))") y_noncrit(1:10)
        write(*,"(*(i0,:,x))") y_crit(1:10)
        write(*,"(*(i0,:,x))") X_noncrit(1,2:11)
        write(*,*) count(X_noncrit(1:10,:)==-1,dim=2)
        call shuffle(X_noncrit, y_noncrit)



    allocate(X_test(N_test, p), X_train(num_noncrit-N_test,p))

        X_test(:,1) = 1
        x_train(:,1) = 1

        X_test(:,2:) = X_noncrit(1:N_test,:)
        X_train(:,2:) = X_noncrit(N_test+1:,:)
        deallocate(X_noncrit)
        y_test = y_noncrit(1:N_test)
        y_train = y_noncrit(N_test+1:)


    simulations: block
        class(binary_logreg), allocatable :: fitter

        real(dp), allocatable :: lambdas(:), momentums(:), learning_rates(:), &
                                 training_accuracies(:,:,:), &
                                 test_accuracies(:,:,:), &
                                 crit_accuracies(:,:,:)
        integer :: num_lambdas, num_momentums, num_learning_rates, i, j, k, iteration
        integer, allocatable :: training_pred(:), test_pred(:), crit_pred(:)

        lambdas = [(10.0d0**i, i = -4,0)]
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
