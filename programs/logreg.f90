program logreg
    use mod_ising
    use mod_utilities, only: shuffle
    use mod_binary_logreg
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, allocatable :: X_crit(:,:), X_noncrit(:,:)
    integer, allocatable :: y_crit(:), y_noncrit(:)
    real(dp), allocatable :: X_train(:,:), X_test(:,:)
    integer, allocatable :: y_test(:), y_train(:), y_pred(:)
    integer :: num_spins, num_ordered, num_crit, num_disordered, num_noncrit
    integer :: u, p

    open(newunit=u, file="data/states.bin", action="read", access="stream")
    read(u) num_spins, num_ordered, num_crit, num_disordered
    write(*,*) num_spins, num_ordered, num_crit, num_disordered

    p = num_spins+1
    num_noncrit = num_ordered + num_disordered

    allocate(X_crit(num_spins,num_crit), y_crit(num_crit), &
             X_noncrit(num_spins,num_noncrit), y_noncrit(num_noncrit))

    read(u) X_noncrit(:,:num_ordered)
    read(u) X_crit(:,:)
    read(u) X_noncrit(:,num_ordered+1:)
    close(u)

    open(newunit=u, file="data/labels.bin", access="stream", status="old")
    read(u) y_noncrit(1:num_ordered), y_crit(:), y_noncrit(num_ordered+1:)
    close(u)

    X_noncrit = transpose(X_noncrit)
    X_crit = transpose(X_crit)
    write(*,*) shape(X_noncrit), shape(y_noncrit)
    write(*,*) shape(X_crit), shape(y_crit)

    write(*,"(*(i0,:,x))") y_noncrit(1:10)
    write(*,"(*(i0,:,x))") y_crit(1:10)
    write(*,"(*(i0,:,x))") X_noncrit(1,2:11)
    write(*,*) count(X_noncrit(1:10,:)==-1,dim=2)
    call shuffle(X_noncrit, y_noncrit)

    setup: block
        real(dp) :: test_fraction
        integer :: N_test

        test_fraction = 0.2d0
        N_test = nint(test_fraction*num_noncrit)

        allocate(X_test(N_test, p), X_train(num_noncrit-N_test,p))

        X_test(:,1) = 1
        x_train(:,1) = 1

        X_test(:,2:) = X_noncrit(1:N_test,:)
        X_train(:,2:) = X_noncrit(N_test+1:,:)
        deallocate(X_noncrit)
        y_test = y_noncrit(1:N_test)
        y_train = y_noncrit(N_test+1:)
    end block setup

    analysis: block
        class(binary_logreg), allocatable :: fitter

        fitter = binary_logreg(learning_rate=0.001d0, &
                               batch_size=200, max_iterations=100, lambda=0.01d0)

        allocate(y_pred, mold=y_test)
        call fitter%fit(X_train, y_train)
        call fitter%predict(X_test, y_pred)

        write(*,*) "Test accuracy:", count(y_pred==y_test)/(1.0d0*size(y_pred))
    end block analysis

end program
