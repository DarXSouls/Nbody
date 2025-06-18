subroutine initialise_run

!
! pcc 18/10/2020
!
! sets lots of parameters before we start
!

use global_variables
use global_parameters
use timing_info
use accel_loop_list

implicit none

integer :: i, nn
integer :: ic_filename_length
character*4 :: filename_tail
integer :: t_power
double precision :: dt_initial(1:N)
integer :: idummy

!
! set the loop controller for the accel function
!
do i = 1, N
   indices(i) = i
end do
do i = N+1, 2*N
   indices(i) = i - N
end do
print *, "INITIALISE: ", indices

!
! set the snapshot time stamps and filecounter
!
time_at_last_snapshot = global_time
time_at_next_snapshot = global_time + time_between_snapshots
if (restart_flag .eq. 1) then
   !
   ! restarted from a snapshot, so need to get the filestamp from the tail
   ! of the ic filename
   !
   ic_filename_length = len_trim( initial_conditions )
   filename_tail = initial_conditions(ic_filename_length-4:ic_filename_length)
   read(filename_tail, *) snap_file_counter
   print *, "Started from a snapshot ending in number ",  snap_file_counter
else
   !
   ! started from a true IC file, so can start the snapshot numbering from 0
   !
   snap_file_counter = 0
end if

!
! Important one... Set G! :) 
!
units_mass_g = solar_mass
units_length_cm = pc
units_time_s = year

G_int = gg * units_mass_g * units_time_s**2 / units_length_cm**3
print *, 'INITIALISE: value of G internally is ', G_int

!
! set the time arrays
!
t(:) = global_time

!
! set the minimum timestep allowed by the integration scheme, 
! based on moving distance soft*eta_dt in dt via interaction with most massive body
!
dt_min = 0.0
!dt_min = sqrt(2.0d0 * 0.1 * soft**3 / (G_int * maxval(m)))
print *, "INITIALISE: forcing minimum dt of ", dt_min
print *, "2.0d0 * eta_dt * soft**3", 2.0d0 * eta_dt * soft**3
print *, "maxval(m)", maxval(m)

!
! get the hermite block time steps scheme up and running, if it's being used.
!
if (int_mode .eq. 2) then

   !
   ! set all particles to be active XXX: this could be removed!
   !
   N_act = N
   do i = 1, N
      i_active(i) = i
   end do

   !
   ! get the accels at the current time
   !
   call get_accel_for_active(N, N_act, i_active, x, y, z, vx, vy, vz, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0, radii)
   a = sqrt(ax0*ax0 + ay0*ay0 + az0*az0)
   adot = sqrt(axdot0*axdot0 + aydot0*aydot0 + azdot0*azdot0)
   dt_initial = min(dt_max, eta_dt * a/adot)
   if ( system_dt .gt. 0 ) then
      do i = 1, N
         if ( is_system(i) .gt. 0)  dt_initial(i) = system_dt
      end do
   end if 

   !
   ! prevent going past output time...
   !
   do i = 1, N 
      if ( (global_time + dt_initial(i)) .gt. time_at_next_snapshot ) then
         dt_initial(i) = time_at_next_snapshot - global_time
      end if
   end do

   !
   ! define log2_dtmax
   !
   log2_dtmax = log(dt_max)/log2

   !
   ! Now convert all dts into a power of 2... 
   !
   t_power_bin = int(log2_dtmax - log(dt_initial)/log2) + 1
   delta_t = dt_max / 2**t_power_bin
   time_bin_pop(:) = 0
   do i = 1, N
      if (t_power_bin(i) .gt. 32) then
         print *, "REALLY SMALL DT ON STARTUP! ABORTING!"
         stop
      end if
      time_bin_pop(t_power_bin(i)) = time_bin_pop(t_power_bin(i)) + 1
   end do
   print *, "Initial delta_t", delta_t
   print *, "Initial t_power_bin", t_power_bin
   print *, "timebin populations: ", time_bin_pop  
 
   ! make the time between snapshots an integer multiple of dt_max
   ! XXX put a print statement here to notify user and check that this works!  
   !t_power = int(log2_dtmax - log(time_between_snapshots)/log2)
   !print *, "INITIALIZE: modifying snap cadence from ", time_between_snapshots
   !print *, 'to ', dt_max / 2**t_power
   !time_between_snapshots = dt_max / 2**t_power
end if

!
! fixed dt!
!
if ( int_mode.eq.0 ) then
   dt = dt_max
   print *, 'FIXED dt INTEGRATION! dt fixed at:', dt
end if

!
! set the timers
!
call system_clock(count_rate=timing_cr)
call system_clock(count_max=timing_cm)
timing_rate = dble(timing_cr)
print *, 'INITIALISE: code timing sample rate in s: ', timing_rate
call system_clock(timing_start)
timing_accel = 0
timing_dts = 0
timing_hermite = 0
timing_herm2 = 0
timing_snaps = 0
timing_diag = 0
accel_function_counter = 0
calls_to_get_accel = 0

end subroutine initialise_run
