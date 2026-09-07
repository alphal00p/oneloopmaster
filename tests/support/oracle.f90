! Independent test driver for the unchanged original OneLOop Fortran library.
! Build (the module and archive paths refer to an existing original build):
! gfortran -I /path/to/OneLOop tests/support/oracle.f90 \
!   /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-oracle
! Protocol: integral kind and unsquared scale, momenta, squared masses, next/stop.
! Kinds 1,2,3,4 are A0,B0,C0,D0; -2 is dB0/dp^2.
program oneloop_oracle
    use avh_olo
    implicit none
    integer :: integral_kind, status, i, legs, np
    real(kind(1d0)) :: mu
    complex(kind(1d0)) :: p(6), m(4), result(0:2)
    character(16) :: command

    do
        read(*, *, iostat=status) integral_kind, mu
        if (status /= 0) exit
        legs = abs(integral_kind)
        select case (integral_kind)
        case (1)
            np = 0
        case (2, -2)
            np = 1
        case (3)
            np = 3
        case (4)
            np = 6
        case default
            error stop 'Unsupported oracle kind'
        end select
        do i = 1, np
            read(*, *) p(i)
        end do
        do i = 1, legs
            read(*, *) m(i)
        end do
        call olo_scale(mu)
        select case (integral_kind)
        case (1)
            call olo(result, m(1))
        case (2)
            call olo(result, p(1), m(1), m(2))
        case (-2)
            call olo_db0(result, p(1), m(1), m(2))
        case (3)
            call olo(result, p(1), p(2), p(3), m(1), m(2), m(3))
        case (4)
            call olo(result, p(1), p(2), p(3), p(4), p(5), p(6), m(1), m(2), m(3), m(4))
        end select
        do i = 0, 2
            write(*, '(a,2es26.17)') 'olo: ', real(result(i)), aimag(result(i))
        end do
        read(*, *, iostat=status) command
        if (status /= 0 .or. trim(command) == 'stop') exit
        if (trim(command) /= 'next') error stop 'Expected next or stop'
    end do
end program oneloop_oracle
