program nbd

!
!
! N-Body Dynamics (NBD) code with at least 4th order Hermite.
!
! Takes an initial condition file and parameter file as input
!

use global_parameters
use global_variables

implicit none

!
!
!
print *, 'precision', tiny(test_quad)

!
! load parameters file
!
call read_param_file

!
! load the initial conditions from file
!
call get_ics

!
! initialise the run
!
call initialise_run

!
! initial diagnostics
!
call diagnostics(0)

!
! The main integration loop
!
integration_loop_counter = 0
do while (global_time .lt. end_time)
   !print *, 'in loop ', integration_loop_counter, global_time, dt 
   !
   ! Is it time for a snapshot? Do energy checks here
   !
   call write_a_snapshot

   !
   ! Do the integration step here. This subroutine takes care of timesteps and
   ! performs calls to get_accel_info.
   !
   if ( int_mode .eq. 0 ) then
      call hermite_4thorder_PECn_fixed_dt
   else if ( int_mode .eq. 1 ) then
      call hermite_4thorder_PECn_global_dt
   else if ( int_mode .eq. 2 ) then
      call hermite_4thorder_PECn_block_dt
   else if ( int_mode .eq. 3 ) then
      call verlet_2ndorder_global_dt
   else
      print *, "invalid parameter int_mode! Must be 0 (fixed), 1 (global), 2 (block), or 3 (verlet, global)"
      stop
   end if

   !
   ! diagnostics
   !
   call diagnostics(1)

   !
   ! update the integration loop counter
   !
   integration_loop_counter = integration_loop_counter + 1
   !if ( integration_loop_counter.eq. 40) stop
end do

end program nbd
