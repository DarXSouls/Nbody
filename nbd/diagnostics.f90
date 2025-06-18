subroutine diagnostics(diagnostic_mode)

!
!
! This code prints out things that the user might like to know, such as information 
! on timesteps, energy, momentum, (V)COM, masses. We've got two modes:
! 	diagnostic_mode = 0 --> initial information from startup / ICs
!	diagnostic_mode = 1 --> information from inside the time integration loop
!

use global_variables
use global_parameters
use timing_info

implicit none

integer :: diagnostic_mode
integer :: i

call get_energies_and_momentum

if ( diagnostic_mode.eq.0 ) then
   !
   ! set the counter to zero
   !
   diagnostic_counter = 0

   !
   ! print out the initial information. Can add to this as more
   ! useful quantities become apparent!
   !
   print *, '------------------------------------------- '
   print *, ' >>>>> the N Body Dynamics (NBD) code <<<<<' 
   print *, '------------------------------------------- '
   print *, 'Using PEC-n scheme with n = ', n_order
   if ( int_mode .eq. 1 ) then
      print *, 'Using global (but variable) timesteps'
   else if (int_mode .eq. 2 ) then
      print *, 'Using block (i.e. variable and individual) timesteps' 
   else
      print *, 'Using a constant timestep for all particles'
   end if 
   print *, ' '
   print *, '>>> Information on INITIAL conditions <<<'
   print *, 'Number of bodies (N): ', N
   print *, 'mass units (g):      ', units_mass_g
   print *, 'distance units (cm): ', units_length_cm
   print *, 'time units (s)     : ', units_time_s
   print *, 'Gravtiational constant: ', G_int
   print *, 'centre of mass:      ', xcom, ycom, zcom
   print *, 'centre of velocity   ', vxcom, vycom, vzcom
   print *, 'total spec angular momentum: ', spec_ang_mom_tot
   print *, 'gravitational and kinetic energies: ', e_grav, e_kin
   print *, 'total energy:                       ', e_grav + e_kin
else if (diagnostic_mode.eq.1) then
   !
   ! update the counter 
   !   
   diagnostic_counter = diagnostic_counter + 1
   !
   ! is it time to write out information?
   !
   if (diagnostic_counter.lt.diagnostic_freq) then
      return
   else
      diagnostic_counter = 0
   end if

   call system_clock(timing_end)
   time_so_far = (timing_end - timing_start)/timing_rate

   !
   ! print out the diagnostics from the loop
   !
   print *, ' '
   print *, 'INFO AT STEP ', integration_loop_counter
   print *, 'time and dt', global_time, dt
   print *, 'centre of mass:      ', xcom, ycom, zcom
   print *, 'centre of velocity   ', vxcom, vycom, vzcom
   print *, 'total spec angular momentum: ', spec_ang_mom_tot
   print *, 'gravitational and kinetic energies: ', e_grav, e_kin
   print *, 'total energy:                       ', e_grav + e_kin      
   print *, 'PERFORMANCE DIAGNOSTICS'
   print *, 'Number of accel pair calculations and calls to get_accel: ', accel_function_counter, calls_to_get_accel 
   print *, 'Total time so far', time_so_far
   print *, 'code time per s (i.e. years/s)', global_time/time_so_far
   if ( int_mode .eq. 2 ) then
      time_bin_pop(:) = 0
      do i = 1, N
         if ( t_power_bin(i) .gt. 32 ) then
            print *, "REALL SMALL DT! ABORTING..."
            stop
         end if 
         time_bin_pop(t_power_bin(i)) = time_bin_pop(t_power_bin(i)) + 1
      end do
      print *, ">> Block timestep population << "
      print *, "Bin number, No. of bodies"
      do i = 1, 32
         print *, i, time_bin_pop(i)
      end do
   end if
   print *, ''
else
   !
   ! doesn't exist!
   !
   print *, "wrong diagnostic mode!"
   stop
end if

end subroutine diagnostics
