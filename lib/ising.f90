module mod_ising
    use mod_basis
    use iso_fortran_env, only: dp => real64
    implicit none

    type, public, extends(basisfunction) :: coupling
        integer :: i, j, L

        contains
            procedure :: eval => eval_coupling
    end type

    type, public, extends(basisfunction) :: spin_function
        integer :: i
        contains
            procedure :: eval => eval_spin
    end type

    contains
        function eval_coupling(self, x) result(H)
            class(coupling), intent(in) :: self
            real(dp), intent(in) :: x(:)
                ! Vectorised matrix of couplings, - si * sj
            real(dp) :: H

            H = -x((self%j-1)*self%L + self%i)
        end function

        function eval_spin(self, x) result(s)
            class(spin_function), intent(in) :: self
            real(dp), intent(in) :: x(:)
            real(dp) :: s

            if (self%i == 0) then
                s = 1.0d0
            else
                s = x(self%i)
            end if
        end function

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

        subroutine create_ising_basis(basis, L)
            class(coupling), allocatable, intent(out) :: basis(:)
            integer, intent(in) :: L

            integer :: i, j

            allocate(basis(L**2))
            do j = 1, L
                do i = 1, L
                    associate(s => basis((j-1)*L + i))
                        s%i = i
                        s%j = j
                        s%L = L
                    end associate
                end do
            end do
        end subroutine

        subroutine create_spin_basis(basis, L)
            class(spin_function), allocatable, intent(out) :: basis(:)
            integer, intent(in) :: L

            integer :: i
            allocate(basis(1:L+1))
            do i = 0, L
                basis(i+1)%i = i
            end do
        end subroutine
end module
