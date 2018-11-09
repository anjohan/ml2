module mod_ising
    use mod_utilities, only: shuffle
    use iso_fortran_env, only: dp => real64
    implicit none

    contains
        function energy(s)
            real(dp), intent(in) :: s(:)
            real(dp) :: energy
            integer :: L

            L = size(s)
            energy = - (s(1)*s(L) + dot_product(s(1:L-1), s(2:L)))
        end function

        function spin_coupling_vector(s) result(x)
            real(dp), intent(in) :: s(:)
            real(dp), allocatable :: x(:)

            integer :: i, j, L
            L = size(s)

            allocate(x(L**2))
            do j = 1, L
                do i = 1, L
                    x((j-1)*L + i) = s(i) * s(j)
                end do
            end do
        end function

        subroutine read_2d_states(states_name, labels_name, &
                                  X_train, X_test, X_crit, &
                                  y_train, y_test, y_crit, &
                                  test_fraction, add_intercept, T)
            character(len=*), intent(in) :: states_name, labels_name
            real(dp), allocatable, intent(out) :: X_train(:,:), X_test(:,:), X_crit(:,:)
            integer, allocatable, intent(out) :: y_train(:), y_test(:), y_crit(:)
            real(dp), intent(in) :: test_fraction
            real(dp), allocatable :: X_noncrit(:,:)
            integer, allocatable :: y_noncrit(:)
            logical, optional, intent(in) :: T
            logical, optional, intent(in) :: add_intercept

            integer :: u, num_spins, num_ordered, num_crit, num_disordered, num_noncrit
            integer :: p, N_test

            open(newunit=u, file=states_name, action="read", access="stream")
            read(u) num_spins, num_ordered, num_crit, num_disordered
            write(*,*) num_spins, num_ordered, num_crit, num_disordered

            num_noncrit = num_ordered + num_disordered
            p = num_spins
            if (present(add_intercept)) p = p + 1

            allocate(X_crit(p,num_crit), &
                X_noncrit(p,num_noncrit), y_noncrit(num_noncrit))

            allocate(y_crit(num_crit))

            if (present(add_intercept)) then
                X_crit(1,:) = 1
                X_noncrit(1,:) = 1
            end if

            read(u) X_noncrit(p-num_spins+1:,:num_ordered)
            read(u) X_crit(p-num_spins+1:,:)
            read(u) X_noncrit(p-num_spins+1:,num_ordered+1:)
            close(u)

            open(newunit=u, file=labels_name, access="stream", status="old")
            read(u) y_noncrit(1:num_ordered), y_crit(:), y_noncrit(num_ordered+1:)
            close(u)

            call shuffle(X_noncrit, y_noncrit, column_wise=.true.)
            call shuffle(X_crit, y_crit, column_wise=.true.)

            N_test = nint(test_fraction*num_noncrit)

            X_test = X_noncrit(:,1:N_test)
            X_train = X_noncrit(:,N_test+1:)
            y_test = y_noncrit(1:N_test)
            y_train = y_noncrit(N_test+1:)
            deallocate(X_noncrit, y_noncrit)

            if (present(T)) then
                X_train = transpose(X_train)
                X_test = transpose(X_test)
                X_crit = transpose(X_crit)
            end if
        end subroutine
end module
