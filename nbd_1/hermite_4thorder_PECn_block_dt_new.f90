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

!
! set the new global_time (i.e. the time at the END of the upcoming step) 
! and various 'dt' arrays
!

call system_clock(th0)

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
   call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1, radii)

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
   do nn = 1, N_act
      i = i_active(nn)
      xp(i) = xp(i)  + t_diff2(i)*axc2(i)/24.0d0 + t_diff2(i)*axc3(i)/120.0d0
      yp(i) = yp(i)  + t_diff2(i)*ayc2(i)/24.0d0 + t_diff2(i)*ayc3(i)/120.0d0
      zp(i) = zp(i)  + t_diff2(i)*azc2(i)/24.0d0 + t_diff2(i)*azc3(i)/120.0d0

      vxp(i) = vxp(i)  + t_diff(i)*axc2(i)/6.0d0 + t_diff(i)*axc3(i)/24.0d0
      vyp(i) = vyp(i)  + t_diff(i)*ayc2(i)/6.0d0 + t_diff(i)*ayc3(i)/24.0d0
      vzp(i) = vzp(i)  + t_diff(i)*azc2(i)/6.0d0 + t_diff(i)*azc3(i)/24.0d0
   end do
end do

!
! update pos / vel for *active* particles with new iterative solution 
!
do nn = 1, N_act
   call system_clock(t0)
   i = i_active(nn)
   x(i) = xp(i)
   y(i) = yp(i)
   z(i) = zp(i)
   vx(i) = vxp(i)
   vy(i) = vyp(i)
   vz(i) = vzp(i)
   !
   ! divide ac2 and ac3 by t for use in the timestep criterion
   !
   recip_t_diff2 = 1.0d0 / t_diff2(i)
   recip_t_diff3 = 1.0d0 / t_diff3(i)
   axc2(i) = axc2(i) * recip_t_diff2
   ayc2(i) = ayc2(i) * recip_t_diff2
   azc2(i) = azc2(i) * recip_t_diff2
   axc3(i) = axc3(i) * recip_t_diff3 
   ayc3(i) = ayc3(i) * recip_t_diff3
   azc3(i) = azc3(i) * recip_t_diff3
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
   ! flip the isok_change_dt flag between pos/neg
   !
   isok_change_dt(i) = isok_change_dt(i) * (-1) 
   call system_clock(t1)
   timing_herm2 = timing_herm2 + (t1 - t0)/timing_rate
 
   !
   ! use the 'mystical' (according to W. Denon!) Aarseth criterion
   ! to set a new dt, if allowed! 
   !
   call system_clock(t0)
   if ( isok_change_dt(i).gt.0 ) then

      a1 = sqrt( ax1(i)*ax1(i) + ay1(i)*ay1(i) + az1(i)*az1(i) )
      a1dot = sqrt( axdot1(i)*axdot1(i) + aydot1(i)*aydot1(i) + azdot1(i)*azdot1(i) ) 
      !ac3 = sqrt( axc3(i)*axc3(i) + ayc3(i)*ayc3(i) + azc3(i)*azc3(i) )
      !acx2_ext = axc2(i) + t_diff(i)*axc3(i)
      !acy2_ext = ayc2(i) + t_diff(i)*ayc3(i)
      !acz2_ext = azc2(i) + t_diff(i)*azc3(i)
      !ac2_ext = sqrt(acx2_ext*acx2_ext + acy2_ext*acy2_ext + acz2_ext*acz2_ext)
      !dt = min(dt_max, sqrt( eta_dt * (a1*ac2_ext + a1dot*a1dot) / (a1dot*ac3 + ac2_ext*ac2_ext) ))
      dt = min(dt_max, eta_dt * a1/a1dot)
      dt = max(dt, 1e-5) 
      !dt = 0.001
      !
      ! Now convert this into a power of 2...
      !
      dt_bin_try = int(log2_dtmax - log(dt)/log2) + 1 
      if (dt_bin_try .gt. t_power_bin(i)) then
         t_power_bin(i) = t_power_bin(i) + 1
      else if ( dt_bin_try .gt. t_power_bin(i) ) then
         t_power_bin(i) = max(0, t_power_bin(i) - 1)
      end if
      delta_t(i) = dt_max / 2**t_power_bin(i)
   end if
   call system_clock(t1)
   timing_dts = timing_dts + (t1 - t0)/timing_rate
end do

!
! get new accel info for the active particles at their new position / velocity
! Note we use the pos / vel values stored in the 'predicted' arrays, since
! these are 'drifted' in space between active times.
!

call get_accel_for_active(N, N_act, i_active, xp, yp, zp, vxp, vyp, vzp, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0, radii)

call system_clock(th1)
timing_hermite = timing_hermite + (th1 - th0)/timing_rate

end subroutine hermite_4thorder_PECn_block_dt
