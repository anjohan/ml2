module mod_functions
    use iso_fortran_env, only: dp => real64
    contains
        subroutine read_states_labels(fname, states, labels)
            character(len=*), intent(in) :: fname
            real(dp), allocatable, intent(out) :: states(:,:), labels(:)
            integer :: num_states, num_spins, u

            open(newunit=u, access="stream", status="old", file="data/" // fname // "_states.bin")
            read(u) num_states, num_spins
            allocate(states(num_states, num_spins))
            read(u) states
            close(u)

            open(newunit=u, access="stream", status="old", file="data/" // fname // "_labels.bin")
            read(u) num_states
            allocate(labels(num_states))
            read(u) labels
            close(u)
        end subroutine
end module

program logreg
    use mod_ising
    use mod_utilities, only: shuffle
    use mod_functions
    use mod_binary_logreg
    use iso_fortran_env, only: dp => real64

    class(binary_logreg), allocatable :: fitter
    class(spin_function), allocatable :: basis(:)

    real(dp), allocatable :: states_disordered(:,:), states_ordered(:,:), states_critical(:,:), &
                             labels_disordered(:), labels_ordered(:), labels_critical(:), &
                             states(:,:), labels(:)

    integer :: num_spins
    integer :: num_disordered, num_ordered, num_critical

    call read_states_labels("ordered", states_ordered, labels_ordered)
    call read_states_labels("disordered", states_disordered, labels_disordered)
    call read_states_labels("critical", states_critical, labels_critical)

    num_spins = size(states_ordered, 2)
    num_ordered = size(states_ordered, 1)
    num_disordered = size(states_disordered, 1)
    num_critical = size(states_critical, 1)

    write(*,*) "Number of spins:", num_spins
    write(*,*) "Shape of ordered:", shape(states_ordered)
    write(*,*) "Shape of disordered:", shape(states_disordered)
    write(*,*) "Shape of critical:", shape(states_critical)

    call shuffle(states_critical, labels_critical)

    write(*,"(*(f0.1,:,x))") labels_ordered(1:10)
    write(*,"(*(f0.1,:,x))") labels_disordered(1:10)
    write(*,"(*(f0.1,:,x))") labels_critical(1:10)


end program
