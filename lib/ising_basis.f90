module mod_ising_basis
    use mod_basis
    implicit none

    type, public, extends(basisfunction) :: coupling
        integer :: i, j, L

        contains
            procedure :: eval => eval_coupling
    end type

    contains
        function eval_coupling(self, x) result(H)
            class(coupling), intent(in) :: self
            real(dp), intent(in) :: x(:)
                ! Vectorised matrix of couplings, - si * sj
            real(dp) :: H

            H = -x((self%j-1)*self%L + self%i)
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
end module
