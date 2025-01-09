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

    subroutine s_mpi_loadbalance_init()

        real(kind(0d0)) :: time_avg
        real(kind(0d0)) :: proc_time(0:num_procs -1)
        logical :: file_exists

        integer :: i, istat, ierr

        if (proc_rank == 0) then
            ! write out the proc_io time data for all the processors
            inquire (FILE='proc_time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='proc_time_data.dat', status='old', action='read')
                read (1, *)
                read(1, '(15F15.8)') (proc_time(i), i = 0, num_procs-1)
                close (1)
            else
                time_avg = -1.0d0
            endif
        end if

        call MPI_BCAST(proc_time, num_procs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        time_avg = proc_time(proc_rank)

        ! if (proc_rank == 0) then
        !     inquire (FILE='repartitioning.dat', EXIST=file_exists)
        !     if (file_exists) then
        !         open(1, file='repartitioning.dat', status='old', action='read')
        !         read(1, *) 
        !         read(1, '(15I5)') proc_coords_x(1:num_procs)
        !         read(1, '(15I5)') proc_coords_y(1:num_procs)
        !         read(1, '(15I5)') proc_coords_z(1:num_procs)
        !         read(1, *) 
        !         read(1, *) 
        !         read(1, '(I5, 9I10)') proc_counts_x 
        !         if (n > 0) then
        !             read(1, '(I5, 9I10)') proc_counts_y
        !             if (p > 0) then
        !                 read(1, '(I5, 9I10)') proc_counts_z
        !             end if
        !         end if
        !         read(1, *)             
        !         close(1)
        !     end if
        ! endif

        if (time_avg > 0.0d0) call s_mpi_loadbalance_computational_domain(time_avg, istat, 40)

    end subroutine s_mpi_loadbalance_init

    !>  The goal of this subroutine is to rebalance the number of points in each direction among the processors
        !! based on the load factor such that the computational load is balanced among the processors.
    subroutine s_mpi_loadbalance_computational_domain(time_avg, istat, opt)
        real(kind(0d0)), intent(inout) :: time_avg
        integer, intent(in), optional :: opt
        integer, intent(out) :: istat

        integer :: i, j, k, ierr
        integer :: tmp_val, mx
        real(kind(0d0)) :: proc_time(0:num_procs -1)
        real(kind(0d0)), dimension(num_procs_x) :: load_factor_x
        real(kind(0d0)), dimension(num_procs_y) :: load_factor_y
        real(kind(0d0)), dimension(num_procs_z) :: load_factor_z
        integer, dimension(num_procs_x) :: px
        integer, dimension(num_procs_y) :: py
        integer, dimension(num_procs_z) :: pz
        integer, dimension(num_procs_x) :: new_dist_x, new_displ_x
        integer, dimension(num_procs_y) :: new_dist_y, new_displ_y
        integer, dimension(num_procs_z) :: new_dist_z, new_displ_z
        integer, dimension(num_procs) :: buffer
        real(kind(0d0)) :: load_factor(num_procs)
        real(kind(0d0)) :: max_load_imbalance
        logical :: file_exists
        integer :: buff_min
        integer, parameter  :: buff_min_threshold = 10

        istat = 1

        buff_min = minval(buff_size_lb)

        ! get the minimum buff_min across all processes 
        if (proc_rank == 0) then
            call MPI_ALLREDUCE(MPI_IN_PLACE, buff_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
        else
            call MPI_ALLREDUCE(buff_min, buff_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
        end if

        if (buff_min < buff_min_threshold .and. .not.present(opt)) then
            return
        end if
        call mpi_bcast_time_step_values(proc_time, time_avg)

        if (proc_rank == 0) then
            ! compute the load imbalance factor
            ! write out the proc_io time data for all the processors
            inquire (FILE='proc_time_data.dat', EXIST=file_exists)
            if (file_exists) then
                open (1, file='proc_time_data.dat', position='append', status='old')
            else
                open (1, file='proc_time_data.dat', status='new')
                write (1, '(A10, A15)') "Ranks", "s/step"
            end if
            write (1, '(15F15.8)') (proc_time(i), i = 0, num_procs-1)
            close (1)
        end if

        diff_start_idx = 0
        diff_count_idx = 0
        mx = buff_size
        if (present(opt)) mx = opt
        ! get the std deviation of the proc_time array
        max_load_imbalance = 0.0d0
        if (proc_rank == 0) then
            load_factor = proc_time / minval(proc_time) - 1.0d0
            max_load_imbalance = maxval(abs(load_factor))
            ! proc_time_std = sqrt(sum((proc_time - sum(proc_time)/num_procs)**2)/num_procs)
            ! proc_time_std = proc_time_std/sum(proc_time)/num_procs
        end if

        call MPI_BCAST(max_load_imbalance, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        ! if the max under utilization is not less than 95 % then return
        if (max_load_imbalance <= 0.0526d0) then
            return
        end if

        if (proc_rank == 0) then
            print *, 'Repartitioning the computational domain'

            load_factor_x = 0.0d0
            load_factor_y = 0.0d0
            load_factor_z = 0.0d0
            do i = 1, num_procs
                load_factor(i) = proc_time(i-1) / minval(proc_time)
                ! X -direction
                do j = 1, num_procs_x
                    if (proc_coords_x(i) == j - 1) then
                        load_factor_x(j) = load_factor_x(j) + load_factor(i)
                        px(j) = proc_counts_x(i) + 1
                    end if
                end do
                if (n == 0) cycle
                ! y -direction
                do j = 1, num_procs_y
                    if (proc_coords_y(i) == j - 1) then
                        py(j) = proc_counts_y(i) + 1
                        load_factor_y(j) = load_factor_y(j) + load_factor(i)
                    end if
                end do
                if (p == 0) cycle
                ! z -direction
                do j = 1, num_procs_z
                    if (proc_coords_z(i) == j - 1) then
                        pz(j) = proc_counts_z(i) + 1
                        load_factor_z(j) = load_factor_z(j) + load_factor(i)
                    end if
                end do
            end do
            load_factor_x = load_factor_x/num_procs_y/num_procs_z
            load_factor_y = load_factor_y/num_procs_x/num_procs_z
            load_factor_z = load_factor_z/num_procs_x/num_procs_y

            call s_redistribute(px, load_factor_x, mx, new_dist_x, new_displ_x)
            if (n > 0) then
                call s_redistribute(py, load_factor_y, mx, new_dist_y, new_displ_y)
            end if
            if (p > 0) then
                call s_redistribute(pz, load_factor_z, mx, new_dist_z, new_displ_z)
            end If
            istat = 0
        end if

        call MPI_BCAST(istat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        !> distribute m, n, and p among the processors
        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_dist_x(proc_coords_x(i) + 1) - 1
            end do
        end if
        call s_mpi_scatter(buffer, tmp_val)
        diff_count_idx(1) = tmp_val - m
        m = tmp_val

        if (n > 0) then
            !
            if (proc_rank == 0) then
                do i = 1, num_procs
                    buffer(i) = new_dist_y(proc_coords_y(i) + 1) - 1
                end do
            end if
            call s_mpi_scatter(buffer, tmp_val)
            diff_count_idx(2) = tmp_val - n
            n = tmp_val

            if ( p > 0) then
                !
                if (proc_rank == 0) then
                    do i = 1, num_procs
                        buffer(i) = new_dist_z(proc_coords_z(i) + 1) - 1
                    end do
                end if
                call s_mpi_scatter(buffer, tmp_val)
                diff_count_idx(3) = tmp_val - p
                p = tmp_val
            end if
        end if

        !> Collect the start_idx values from all the processors

        if (proc_rank == 0) then
            do i = 1, num_procs
                buffer(i) = new_displ_x(proc_coords_x(i) + 1) 
            end do
        end if
        call s_mpi_scatter(buffer, tmp_val)
        diff_start_idx(1) = tmp_val - start_idx(1)
        start_idx(1) = tmp_val

        if (n > 0) then
            if (proc_rank == 0) then
                do i = 1, num_procs
                    buffer(i) = new_displ_y(proc_coords_y(i) + 1) 
                end do
            end if
            call s_mpi_scatter(buffer, tmp_val)
            diff_start_idx(2) = tmp_val - start_idx(2)
            start_idx(2) = tmp_val

            if (p > 0) then
                if (proc_rank == 0) then
                    do i = 1, num_procs
                        buffer(i) = new_displ_z(proc_coords_z(i) + 1) 
                    end do
                end if
                call s_mpi_scatter(buffer, tmp_val)
                diff_start_idx(3) = tmp_val - start_idx(3)
                start_idx(3) = tmp_val
            end if
        end if

        call MPI_ALLGATHER(m, 1, MPI_INTEGER, proc_counts_x, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

        if (n > 0) then
            call MPI_ALLGATHER(n, 1, MPI_INTEGER, proc_counts_y, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

            if (p > 0) then
                call MPI_ALLGATHER(p, 1, MPI_INTEGER, proc_counts_z, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
            endif
        endif

        if (proc_rank == 0) then
            inquire (FILE='repartitioning.dat', EXIST=file_exists)
            if (file_exists) then
                open(1, file='repartitioning.dat', status='old', position='append')
            else
                open(1, file='repartitioning.dat', status='new')
                write(1, '(I5)') num_procs
                write(1, '(15I5)') proc_coords_x(1:num_procs)
                write(1, '(15I5)') proc_coords_y(1:num_procs)
                write(1, '(15I5)') proc_coords_z(1:num_procs)
            endif
            write(1, '(A, I5)') 'buff min: ', buff_min
            write(1, '(A, F10.5)') 'max. load imbalance: ', max_load_imbalance
            write(1, '(I5, 9I10)') proc_counts_x 
            if (n > 0) then
                write(1, '(I5, 9I10)') proc_counts_y
                if (p > 0) then
                    write(1, '(I5, 9I10)') proc_counts_z
                end if
            end if
            write(1, *) '----------------------------------------'
            close(1)
        end if

        call s_mpi_barrier()
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
#ifdef DEBUG        
        if (sum(new_dist) /= sum(counts)) then
            print *, 'Error: sum(new_dist) /= sum(counts)'
            print *, 'sum(new_dist) = ', sum(new_dist)
            print *, 'sum(counts) = ', sum(counts)
            error stop 'Error: sum(new_dist) /= sum(counts)'
        end if
#endif        
        deallocate (new_dist_real)
    end subroutine s_redistribute

end module m_load_balance

