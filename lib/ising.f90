module mod_ising
    use mod_basis
    use iso_fortran_env, only: dp => real64
    implicit none

    type, public, extends(basisfunction) :: coupling
        integer :: i, j, L

        contains
            procedure :: eval
    end type

    contains
        function eval(self, x) result(H)
            class(coupling), intent(in) :: self
            real(dp), intent(in) :: x(:)
                ! Vectorised matrix of couplings, si * sj
            real(dp) :: H

            H = x((self%j-1)*self%L + self%i)
        end function

        function energy(s)
            real(dp), intent(in) :: s(:)
            real(dp) :: energy
            integer :: L

            L = size(s)
            energy = s(1)*s(L) + dot_product(s(1:L-1), s(2:L))
        end function
end module
