module m_load_balance

    ! Dependencies =============================================================
#ifdef MFC_MPI
    use mpi                    !< Message passing interface (MPI) module
#endif

    use m_helper_basic         !< Functions to compare floating point numbers

    use m_helper

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_common

    use m_nvtx

    use ieee_arithmetic

    implicit none

contains
    !>  The goal of this subroutine is to rebalance the number of points in each direction among the processors
        !! based on the load factor such that the computational load is balanced among the processors.
    subroutine s_mpi_loadbalance_computational_domain(time_avg)
        real(kind(0d0)), intent(in) :: time_avg
        integer :: i, j, k, ierr, tmp_val

        real(kind(0d0)) :: proc_time(0:num_procs -1)
        real(kind(0d0)), dimension(num_procs_x) :: load_factor_x
        real(kind(0d0)), dimension(num_procs_y) :: load_factor_y
        integer, dimension(num_procs_x) :: px
        integer, dimension(num_procs_y) :: py
        integer, dimension(num_procs_x) :: new_dist_x, new_displ_x
        integer, dimension(num_procs_y) :: new_dist_y, new_displ_y
        integer, dimension(num_procs) :: buffer

        call mpi_bcast_time_step_values(proc_time, time_avg)

        if (time_avg < epsilon(0d0)) then
            return
        end if

        if (proc_rank == 0) then
            load_factor_x = 0.0d0
            load_factor_y = 0.0d0
            do i = 1, num_procs
                load_factor(i) = proc_time(i-1) / minval(proc_time)
                do j = 1, num_procs_x
                    if (proc_coords_x(i) == j - 1) then
                        load_factor_x(j) = load_factor_x(j) + load_factor(i)
                        px(j) = proc_counts_x(i) + 1
                    end if
                end do
                do j = 1, num_procs_y
                    if (proc_coords_y(i) == j - 1) then
                        py(j) = proc_counts_y(i) + 1
                        load_factor_y(j) = load_factor_y(j) + load_factor(i)
                    end if
                end do
            end do
            load_factor_x = load_factor_x/num_procs_y
            load_factor_y = load_factor_y/num_procs_x

            call s_redistribute(px, load_factor_x, 40, new_dist_x, new_displ_x)
            call s_redistribute(py, load_factor_y, 40, new_dist_y, new_displ_y)
        end if
        call s_mpi_barrier()

        !> distribute m and n
        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_dist_x(proc_coords_x(i) + 1) - 1
            end do
        end if
        call s_mpi_scatter(buffer, m)

        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_dist_y(proc_coords_y(i) + 1) - 1
            end do
        end if
        call s_mpi_scatter(buffer, n)

        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_displ_x(proc_coords_x(i) + 1) 
            end do
        end if
        call s_mpi_scatter(buffer, tmp_val)
        diff_start_idx(1) = tmp_val - start_idx(1)
        start_idx(1) = tmp_val

        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_displ_y(proc_coords_y(i) + 1) 
            end do
        end if
        call s_mpi_scatter(buffer, tmp_val)
        diff_start_idx(2) = tmp_val - start_idx(2)
        start_idx(2) = tmp_val

        if (proc_rank == 0) then
            open(1, file='repartitioning.dat', status='new', action='write')
            write(1, '(I5)') num_procs
            do i = 1, num_procs
                write(1, '(3I5)') i, proc_coords_x(i), proc_coords_y(i)
            end do
            close(1)
        end if

    end subroutine s_mpi_loadbalance_computational_domain

    subroutine s_redistribute(counts, lf, mx, new_dist, new_displ)
        implicit none
        integer, intent(in) :: mx
        integer, intent(in) :: counts(:)
        real(kind(0d0)), intent(in) :: lf(:)
        integer, intent(out) :: new_dist(size(counts)), new_displ(size(counts))
        integer :: i, nsz, diff, sum_counts, sum_new_dist
        real(kind(0d0)) :: cost_uniform, scaling_factor
        integer :: mx_allowed

        ! Local variables for calculations
        real(kind=8), allocatable :: new_dist_real(:)
        integer :: max_diff

        nsz = size(counts)
        mx_allowed = mx - 1

        allocate (new_dist_real(nsz))

        ! Compute cost_uniform
        cost_uniform = sum(lf)/real(nsz)

        ! Compute new_dist as a real array
        do i = 1, nsz
            new_dist_real(i) = cost_uniform*real(counts(i))/lf(i)
        end do

        ! Scale new_dist_real to ensure integer totals match counts
        scaling_factor = sum(real(counts))/sum(new_dist_real)
        do i = 1, nsz
            new_dist(i) = floor(new_dist_real(i)*scaling_factor)
        end do

        ! Adjust difference between counts and new_dist
        sum_counts = sum(counts)
        sum_new_dist = sum(new_dist)
        diff = sum_counts - sum_new_dist

        if (diff > 0) then
            do i = 1, nsz
                new_dist(i) = new_dist(i) + 1
            end do
        elseif (diff < 0) then
            do i = 1, nsz
                new_dist(i) = new_dist(i) - 1
            end do
        end if

        ! Ensure that no element deviates by more than mx_allowed from original counts
        do i = 1, nsz
            if (new_dist(i) > counts(i) + mx_allowed) then
                new_dist(i) = counts(i) + mx_allowed
            elseif (new_dist(i) < counts(i) - mx_allowed) then
                new_dist(i) = counts(i) - mx_allowed
            end if
        end do

        ! Recalculate the difference after adjustments
        sum_new_dist = sum(new_dist)
        diff = sum_counts - sum_new_dist

        if (diff > 0) then
            do i = 1, nsz
                if (new_dist(i) < counts(i) + mx_allowed) then
                    new_dist(i) = new_dist(i) + 1
                    diff = diff - 1
                end if
                if (diff == 0) exit
            end do
        elseif (diff < 0) then
            do i = 1, nsz
                if (new_dist(i) > counts(i) - mx_allowed) then
                    new_dist(i) = new_dist(i) - 1
                    diff = diff + 1
                end if
                if (diff == 0) exit
            end do
        end if

        ! Re-scale new_dist to ensure integer totals match counts
        scaling_factor = sum(real(counts))/sum(real(new_dist))
        do i = 1, nsz
            new_dist(i) = floor(real(new_dist(i))*scaling_factor)
        end do

        ! Final adjustment for difference
        diff = sum_counts - sum(new_dist)
        do i = 1, diff
            new_dist(i) = new_dist(i) + 1
        end do

        ! Compute displacements (new_displ)
        new_displ(1) = 0
        do i = 2, nsz
            new_displ(i) = new_displ(i - 1) + new_dist(i - 1)
        end do
        ! Compute maximum difference for assertion
        max_diff = 0
        do i = 1, nsz
            max_diff = max(max_diff, abs(new_dist(i) - counts(i)))
        end do
        deallocate (new_dist_real)
    end subroutine s_redistribute

    subroutine s_balance_elements(N_, size, counts, displ)
        integer, intent(in) :: N_, size
        integer, intent(out) :: counts(:), displ(:)
        integer :: i, count, diff

        ! Calculate base count per division
        count = floor(N_*1.0d0/size)
        counts = count
        diff = N_ - count*size

        ! Adjust the first 'diff' elements of counts
        do i = 1, diff
            counts(i) = counts(i) + 1
        end do

        ! Calculate displacements
        displ(1) = 0
        do i = 2, size
            displ(i) = displ(i - 1) + counts(i - 1)
        end do
    end subroutine s_balance_elements

end module m_load_balance

