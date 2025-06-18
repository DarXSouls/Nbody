subroutine write_a_snapshot
    ! Use necessary modules
    use global_variables
    use global_parameters
    use timing_info

    implicit none
    character*4 :: filename_tail
    integer :: snapfile_len
    integer :: i
    double precision :: dx, dy, dz
    double precision :: dvx, dvy, dvz
    double precision :: v_rel2, r_rel, rdotv
    double precision :: e_bin, a_bin, p_bin, energy_bin
    double precision :: m_sys, m_red
    double precision :: t_predict, snapshot_time

    ! New variables for directory handling
    character(len=100) :: output_directory  
    character(len=100) :: base_directory
    integer :: system_status
    character(len=200) :: command

    
    

    ! Get the current system clock
    call system_clock(tw0)

    ! Generate a unique directory name using the existing tw0 variable
    write(base_directory, '(I10)') tw0
    output_directory = "/home/dxs_cosmos/bin3/output_" // trim(adjustl(base_directory)) 

    ! Create the new directory
    command = "mkdir -p " // trim(output_directory)
    call execute_command_line(trim(command), wait=.true., exitstat=system_status)
    if (system_status /= 0) then
        print *, "Error creating directory: ", trim(output_directory)
        stop
    endif

    ! Find the time that our next integration step will take us to 
    if (int_mode .eq. 2) then
        t_predict = minval(t + delta_t)
    else
        t_predict = global_time + dt
    end if

    ! Decide whether or not we need to write a snapshot
    if ( t_predict .ge. time_at_next_snapshot ) then
        print *, "WRITESNAP, T-PREDICT:", t_predict, time_at_next_snapshot

        ! Set the filename and open
1000    format(I4.4)
        write(filename_tail, 1000) snap_file_counter
        snapfile_len = len_trim(snapshot_base)
        snapshot_filename = trim(output_directory) // "/" // trim(snapshot_base) // filename_tail
        open(30, file=snapshot_filename, form='unformatted')

        ! Synchronize all the particles at the output time
        if (int_mode.eq.0  .or.  int_mode.eq.3) then
            xo(:) = x(:)
            yo(:) = y(:)
            zo(:) = z(:)
            vxo(:) = vx(:)
            vyo(:) = vy(:)
            vzo(:) = vz(:)
            snapshot_time = global_time
        else
            snapshot_time = time_at_next_snapshot
            call synchronise_posvel_for_output
        end if

        call get_energies_and_momentum

        ! Notify user
        print *, 'WRITING SNAPSHOT to file ', snapshot_filename, ' at time ', snapshot_time

        ! Write the data and close
        write(30) N
        write(30) snapshot_time
        write(30) units_mass_g, units_length_cm, units_time_s
        !write(30) G_int
        write(30) xo(1:N)
        write(30) yo(1:N)
        write(30) zo(1:N)
        write(30) vxo(1:N)
        write(30) vyo(1:N)
        write(30) vzo(1:N)
        write(30) m(1:N) 
        write(30) radii(1:N, 1:N)
        close(30)

        ! Write an ASCII file if N is small
        if (N.gt.2  .and.  N.lt.30) then
            if ( snap_file_counter.eq.0 ) then
                open(40, file=trim(output_directory) // "/body_evolution.dat")
            end if
            write(40, *) N
            write(40, *) snapshot_time
            do i = 1, N
                write(40, 4000) xo(i),yo(i),zo(i),vxo(i),vyo(i),vzo(i),m(i)
            end do
        end if

        ! Write out pos, vel, m if N is small
        if (N.le.10) then
            if ( snap_file_counter.eq.0 ) then
                open(50, file=trim(output_directory) // "/orbit_evolution.dat")
                open(60, file=trim(output_directory) // "/binary_parameters.dat")
            end if
            write(50, 4000) snapshot_time,xo(1),yo(1), zo(1), xo(2),yo(2), zo(2)
            dx = xo(2) - xo(1)
            dy = yo(2) - yo(1)
            dz = zo(2) - zo(1)
            dvx = vxo(2) - vxo(1)
            dvy = vyo(2) - vyo(1)
            dvz = vzo(2) - vzo(1)
            r_rel = dsqrt(dx*dx + dy*dy + dz*dz)
            v_rel2 = dvx*dvx + dvy*dvy + dvz*dvz
            m_sys = m(1) + m(2)
            m_red = (m(1) * m(2)) / m_sys
            energy_bin = 0.5d0 * m_red * v_rel2 - G_int * m(1) * m(2) / r_rel
            a_bin = - G_int * m(1) * m(2) / 2.0d0 / energy_bin
            p_bin = dsqrt(a_bin**3 / m_sys)
            rdotv = dx*dvx + dy*dvy + dz*dvz
            e_bin = dsqrt((1.0d0 - r_rel/a_bin)**2 + (rdotv**2)/a_bin/G_int/m_sys)
            write(60, 4000) snapshot_time, dt, energy_bin, a_bin, p_bin, e_bin, r_rel, sqrt(v_rel2)
        end if

4000 format(16(E20.14, 1X))

        ! Update the time stamps and snapcounter
        time_at_last_snapshot = time_at_next_snapshot
        time_at_next_snapshot = time_at_last_snapshot + time_between_snapshots
        print *, 'Next snapshot will be written at time ', time_at_next_snapshot

        ! Update the snapshot file counter
        snap_file_counter = snap_file_counter + 1
    end if

    ! Update timing information
    call system_clock(tw1)
    timing_snaps = timing_snaps + (tw1-tw0)/timing_rate

end subroutine write_a_snapshot
