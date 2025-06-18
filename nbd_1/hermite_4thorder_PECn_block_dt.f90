subroutine hermite_4thorder_PECn_block_dt

use global_variables
use global_parameters
use timing_info

implicit none

double precision :: t_diff(1:N), t_diff2(1:N), t_diff3(1:N)
double precision :: recip_t_diff2, recip_t_diff3
integer :: n_iter
integer :: i, nn
integer :: dt_bin_try
integer :: n_increased_dt, n_decreased_dt
double precision :: delta_t_try

n_increased_dt = 0
n_decreased_dt = 0

!
! set the new global_time (i.e. the time at the END of the upcoming step) 
! and various 'dt' arrays
!

global_time = minval(t + delta_t)
t_diff = global_time - t
t_diff2 = t_diff*t_diff
t_diff3 = t_diff*t_diff*t_diff
dt = minval(t_diff)

!
! Now get a list of active particles (i.e. those that have t = global_time)
!
N_act = 0
do i = 1, N
   if ( abs(t(i) + delta_t(i) - global_time) < 2*tiny(global_time) ) then
      N_act = N_act + 1
      i_active(N_act) = i
   end if
end do
if( N_act .eq. 0) then
   print *, 'PROBLEM: no active particles! Aborting! :-/'
   print *, 'global time, dt', global_time, dt
   print *, 'delta_t', delta_t
   print *, 't_diff', t_diff
   print *, 'abs(t(i) + delta_t(i) - global_time)', abs(t + delta_t - global_time)
   print *, 'time loop counter', integration_loop_counter
   stop
end if

if ( mod(integration_loop_counter, 1000).eq.0 ) then
   print *, "List of active particles at loop", integration_loop_counter
   print *, i_active
   print *, "global time:", global_time
   print *, "active particle times: ", t(i_active(:))
   print *, "active particle dts: ", delta_t(i_active(:))
end if

!
! using the current position and velocity, predict (p) new point
!
xp = x + t_diff*vx + 0.5d0*t_diff2*ax0 + sixth*t_diff3*axdot0
yp = y + t_diff*vy + 0.5d0*t_diff2*ay0 + sixth*t_diff3*aydot0
zp = z + t_diff*vz + 0.5d0*t_diff2*az0 + sixth*t_diff3*azdot0

vxp = vx + t_diff*ax0 + 0.5d0*t_diff2*axdot0
vyp = vy + t_diff*ay0 + 0.5d0*t_diff2*aydot0
vzp = vz + t_diff*az0 + 0.5d0*t_diff2*azdot0

!
! apply the iterative corrector loop as discussed in Kokubo et al. 1998.
! NOTE: this is only really worth it when using a constant global timestep
! but we've kept it here for testing purposes
!
do n_iter = 1, n_order
   !
   ! get the accels based on the predicted or corrected values:
   !
   if (n_iter .eq. 1 ) then
      call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1, radii)
   else
      call get_accel_for_active(N, N_act, i_active, xc, yc, zc, vxc, vyc, vzc, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1, radii)
   end if

   !
   ! now work out the 2nd and 3rd order "corrected" accels, which are based on a0 and a1:
   ! NOTE1: we don't divide by t here, to prevent rounding error in correction)
   !        but rather do it at the end of the step, before working out new dts 
   !
   do nn = 1, N_act
      i = i_active(nn)
      axc2(i) = -6.0d0*(ax0(i) - ax1(i)) - t_diff(i)*(4.0d0*axdot0(i) + 2.0d0*axdot1(i))
      ayc2(i) = -6.0d0*(ay0(i) - ay1(i)) - t_diff(i)*(4.0d0*aydot0(i) + 2.0d0*aydot1(i))
      azc2(i) = -6.0d0*(az0(i) - az1(i)) - t_diff(i)*(4.0d0*azdot0(i) + 2.0d0*azdot1(i))

      axc3(i) = 12.0d0*(ax0(i) - ax1(i)) + 6.0d0*t_diff(i)*(axdot0(i) + axdot1(i))
      ayc3(i) = 12.0d0*(ay0(i) - ay1(i)) + 6.0d0*t_diff(i)*(aydot0(i) + aydot1(i))
      azc3(i) = 12.0d0*(az0(i) - az1(i)) + 6.0d0*t_diff(i)*(azdot0(i) + azdot1(i))
   end do

   !
   ! finally, correct the predicted step using these corrected accels
   !
   xc = xp
   yc = yp
   zc = zp
   vxc = vxp
   vyc = vyp
   vzc = vzp
   do nn = 1, N_act
      i = i_active(nn)
      xc(i) = xp(i)  + t_diff2(i)*axc2(i)/24.0d0 + t_diff2(i)*axc3(i)/120.0d0
      yc(i) = yp(i)  + t_diff2(i)*ayc2(i)/24.0d0 + t_diff2(i)*ayc3(i)/120.0d0
      zc(i) = zp(i)  + t_diff2(i)*azc2(i)/24.0d0 + t_diff2(i)*azc3(i)/120.0d0

      vxc(i) = vxp(i)  + t_diff(i)*axc2(i)/6.0d0 + t_diff(i)*axc3(i)/24.0d0
      vyc(i) = vyp(i)  + t_diff(i)*ayc2(i)/6.0d0 + t_diff(i)*ayc3(i)/24.0d0
      vzc(i) = vzp(i)  + t_diff(i)*azc2(i)/6.0d0 + t_diff(i)*azc3(i)/24.0d0
   end do
end do

!
! update pos / vel for *active* particles with new iterative solution and set the 'predicted' values
! to these too, for use in the accel function
!
do nn = 1, N_act
   i = i_active(nn)
   x(i) = xc(i)
   y(i) = yc(i)
   z(i) = zc(i)
   vx(i) = vxc(i)
   vy(i) = vyc(i)
   vz(i) = vzc(i)
   !
   ! update the drifted / predicted arrays too (for the accel calculation below)
   !
   xp(i) = xc(i)
   yp(i) = yc(i)
   zp(i) = zc(i)
   vxp(i) = vxc(i)
   vyp(i) = vyc(i)
   vzp(i) = vzc(i)
   !
   ! divide ac2 and ac3 by t for use in the timestep criterion
   !
   !recip_t_diff2 = 1.0d0 / t_diff2(i)
   !recip_t_diff3 = 1.0d0 / t_diff3(i)
   !axc2(i) = axc2(i) * recip_t_diff2
   !ayc2(i) = ayc2(i) * recip_t_diff2
   !azc2(i) = azc2(i) * recip_t_diff2
   !axc3(i) = axc3(i) * recip_t_diff3 
   !ayc3(i) = ayc3(i) * recip_t_diff3
   !azc3(i) = azc3(i) * recip_t_diff3
   !
   ! update time and  get new timesteps for ACTIVE particles. Note that we
   ! don't prevent particles evolving past the output time here. Instead, 
   ! we recalculate the positions for output separately.
   ! 
   !
   ! set time of particle to global_time + current dt
   !
   t(i) = t(i) + delta_t(i)

   !
   ! using the simple timestep criterion from Aarseth
   ! to set a new dt, if allowed! 
   !
   if ( is_system(i).eq.0 ) then
      a1 = sqrt( ax1(i)*ax1(i) + ay1(i)*ay1(i) + az1(i)*az1(i) )
      a1dot = sqrt( axdot1(i)*axdot1(i) + aydot1(i)*aydot1(i) + azdot1(i)*azdot1(i) ) 
      dt = min(dt_max, eta_dt * a1/a1dot)
      dt = max(dt, 1e-5) 
      
      !
      ! Now convert this into a power of 2... Integer power can always get larger (smaller dt),
      ! but can only go down after two steps.
      !
      dt_bin_try = int(log2_dtmax - log(dt)/log2) + 1 
      delta_t_try = dt_max / 2**dt_bin_try
      if (dt_bin_try .gt. t_power_bin(i)) then
         t_power_bin(i) = dt_bin_try
         delta_t(i) = delta_t_try
         n_decreased_dt = n_decreased_dt + 1
      else if ( dt_bin_try .lt. t_power_bin(i) ) then
         if ( dmod(t(i) + delta_t_try, 2.0D0).eq.0 ) then
            t_power_bin(i) = max(0, t_power_bin(i) - 1)
            delta_t(i) = dt_max / 2**t_power_bin(i)
            n_increased_dt = n_increased_dt + 1
        end if
      end if
      !print *, "dt_bin_try, delta_t_try, mod term", dt_bin_try, delta_t_try, dmod(t(i) + delta_t_try, 2.0D0)

      !if (dt_bin_try .gt. t_power_bin(i)) then
      !   t_power_bin(i) = dt_bin_try
      !   isok_change_dt(i) = 0
      !else if ( dt_bin_try .lt. t_power_bin(i) .and. isok_change_dt(i).eq.2 ) then
      !   t_power_bin(i) = max(0, t_power_bin(i) - 1)
      !   isok_change_dt(i) = 0
      !end if
      !delta_t(i) = dt_max / 2**t_power_bin(i)
   end if
end do

!
! get new accel info for the active particles at their new position / velocity
! Note we use the pos / vel values stored in the 'predicted' arrays, since
! these are 'drifted' in space between active times.
!

call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0, radii)

!print *, "timesteps up / down this step: ", n_increased_dt, n_decreased_dt

end subroutine hermite_4thorder_PECn_block_dt
