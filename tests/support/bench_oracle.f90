! Optional throughput driver for an unchanged original OneLOop DP build.
! gfortran -O2 -I /path/to/OneLOop tests/support/bench_oracle.f90 &
!   /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-bench
! Input: OLO_BENCH_V1, then nfixtures warmup calls repetitions, followed by
! nfixtures ordinary numeric rows from tests/data/benchmark.txt (no comments).
! SAME repeats one row; HETERO cycles one family's rows in input order.
! Every call computes all three Laurent coefficients. Only the finite complex
! coefficient is accumulated inside timing; full outputs are checked outside.
program oneloop_bench_oracle
    use avh_olo, only: olo, olo_db0
    use iso_fortran_env, only: int64, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    type fixture
        integer :: kind
        real(real64) :: mu
        complex(real64) :: p(6), m(4)
    end type fixture
    type(fixture), allocatable :: rows(:)
    integer, allocatable :: selection(:)
    integer, parameter :: families(5) = [1, 2, -2, 3, 4]
    integer :: nrows, repetitions, row, family, nselected, np, nm, nfields, j, status
    integer(int64) :: warmup, calls, clock_rate, clock_max
    real(real64) :: values(22)
    complex(real64) :: result(0:2)
    character(4096) :: line

    read(*, '(a)', iostat=status) line
    if (status /= 0 .or. trim(line) /= 'OLO_BENCH_V1') error stop 'Expected OLO_BENCH_V1'
    read(*, *, iostat=status) nrows, warmup, calls, repetitions
    if (status /= 0) error stop 'Expected fixture/warmup/call/repetition counts'
    if (nrows < 1 .or. nrows > 100000 .or. warmup < 0 .or. calls < 1 .or. repetitions < 1) &
        error stop 'Invalid benchmark counts'
    allocate(rows(nrows), selection(nrows))
    do row = 1, nrows
        read(*, '(a)', iostat=status) line
        if (status /= 0) error stop 'Missing benchmark fixture'
        read(line, *, iostat=status) rows(row)%kind
        if (status /= 0) error stop 'Missing fixture kind'
        call dimensions(rows(row)%kind, np, nm)
        nfields = 2 + np + 2*nm + 6
        read(line, *, iostat=status) values(1:nfields)
        if (status /= 0) error stop 'Malformed benchmark fixture'
        if (.not. all(ieee_is_finite(values(1:nfields))) .or. values(2) <= 0) &
            error stop 'Nonfinite fixture or invalid squared scale'
        rows(row)%mu = sqrt(values(2))
        rows(row)%p = cmplx(0.0_real64, 0.0_real64, real64)
        rows(row)%m = cmplx(0.0_real64, 0.0_real64, real64)
        do j = 1, np
            rows(row)%p(j) = cmplx(values(2+j), 0.0_real64, real64)
        end do
        do j = 1, nm
            rows(row)%m(j) = cmplx(values(1+np+2*j), values(2+np+2*j), real64)
        end do
        call evaluate(rows(row), result)
        write(*, '(a,1x,i0,1x,i0,6(1x,es26.17))') &
            'VALUE', row, rows(row)%kind, (real(result(j), real64), aimag(result(j)), j=0,2)
    end do
    call system_clock(count_rate=clock_rate, count_max=clock_max)
    if (clock_rate <= 0) error stop 'Unavailable wall clock'
    write(*, '(a,1x,i0,1x,i0)') 'CLOCK', clock_rate, clock_max

    do row = 1, nrows
        selection(1) = row
        call workload('SAME', rows(row)%kind, row, selection(1:1))
    end do
    do family = 1, size(families)
        nselected = 0
        do row = 1, nrows
            if (rows(row)%kind /= families(family)) cycle
            nselected = nselected + 1
            selection(nselected) = row
        end do
        if (nselected > 0) call workload('HETERO', families(family), 0, selection(1:nselected))
    end do

contains

    subroutine dimensions(kind, momenta, masses)
        integer, intent(in) :: kind
        integer, intent(out) :: momenta, masses
        masses = abs(kind)
        select case (kind)
        case (1)
            momenta = 0
        case (2, -2)
            momenta = 1
        case (3)
            momenta = 3
        case (4)
            momenta = 6
        case default
            error stop 'Unsupported scalar benchmark kind'
        end select
    end subroutine dimensions

    subroutine evaluate(input, output)
        type(fixture), intent(in) :: input
        complex(real64), intent(out) :: output(0:2)
        ! The explicit rmu overload receives an UNSQUARED, precomputed scale.
        ! Do not add olo_scale or a square root to the timed call path.
        select case (input%kind)
        case (1)
            call olo(output, input%m(1), input%mu)
        case (2)
            call olo(output, input%p(1), input%m(1), input%m(2), input%mu)
        case (-2)
            call olo_db0(output, input%p(1), input%m(1), input%m(2), input%mu)
        case (3)
            call olo(output, input%p(1), input%p(2), input%p(3), &
                input%m(1), input%m(2), input%m(3), input%mu)
        case (4)
            call olo(output, input%p(1), input%p(2), input%p(3), input%p(4), input%p(5), input%p(6), &
                input%m(1), input%m(2), input%m(3), input%m(4), input%mu)
        end select
    end subroutine evaluate

    subroutine workload(mode, kind, fixture_id, indices)
        character(*), intent(in) :: mode
        integer, intent(in) :: kind, fixture_id, indices(:)
        integer :: cursor, repetition
        integer(int64) :: iteration, start_tick, end_tick
        real(real64) :: wall_seconds, cpu_start, cpu_end
        complex(real64) :: output(0:2), checksum

        cursor = 1
        do iteration = 1, warmup
            call evaluate(rows(indices(cursor)), output)
            cursor = cursor + 1
            if (cursor > size(indices)) cursor = 1
        end do
        do repetition = 1, repetitions
            cursor = 1
            checksum = cmplx(0.0_real64, 0.0_real64, real64)
            call cpu_time(cpu_start)
            call system_clock(start_tick)
            do iteration = 1, calls
                call evaluate(rows(indices(cursor)), output)
                checksum = checksum + output(0)
                cursor = cursor + 1
                if (cursor > size(indices)) cursor = 1
            end do
            call system_clock(end_tick)
            call cpu_time(cpu_end)
            if (end_tick >= start_tick) then
                wall_seconds = real(end_tick-start_tick, real64)/real(clock_rate, real64)
            else
                wall_seconds = (real(clock_max-start_tick, real64)+real(end_tick, real64)+1)/real(clock_rate, real64)
            end if
            write(*, '(a,1x,a,5(1x,i0),10(1x,es26.17))') 'SAMPLE', mode, kind, fixture_id, &
                repetition, calls, size(indices), wall_seconds, cpu_end-cpu_start, &
                real(checksum, real64), aimag(checksum), (real(output(j), real64), aimag(output(j)), j=0,2)
        end do
    end subroutine workload
end program oneloop_bench_oracle
