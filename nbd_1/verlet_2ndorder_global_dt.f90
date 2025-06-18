subroutine verlet_2ndorder_global_dt

use global_variables
use global_parameters

implicit none

double precision :: dt_try(1:N)
integer :: i

!
! get the accels at the current time
!

call get_accel_info(N, x, y, z, vx, vy, vz, m, soft, ax0, ay0, az0, axdot0, aydot0, azdot0, radii)

if ( integration_loop_counter.eq. 0 ) then
   print *, 'Initiallising dts for global time-stepping algortithm' 
   do i = 1, N
      a = sqrt(ax0(i)*ax0(i) + ay0(i)*ay0(i) + az0(i)*az0(i))
      adot = sqrt(axdot0(i)*axdot0(i) + aydot0(i)*aydot0(i) + azdot0(i)*azdot0(i))
      dt_try(i) = eta_dt * a(i)/adot(i)
   end do
   dt = minval(dt_try)
   dt = min(dt, dt_max)
   if ( (global_time + dt) .gt. time_at_next_snapshot ) then
      dt = time_at_next_snapshot - global_time
      print *, "VERLET: dt limited by snapshot time at startup!", global_time, dt
   end if   
end if

!
! using the current position and velocity, get new position
!
x = x + dt*vx + 0.5d0*dt*dt*ax0
y = y + dt*vy + 0.5d0*dt*dt*ay0
z = z + dt*vz + 0.5d0*dt*dt*az0

!
! now use the updated positions to get a new acceleration
!
call get_accel_info(N, x, y, z, vx, vy, vz, m, soft, ax1, ay1, az1, axdot1, aydot1, azdot1, radii)

vx = vx + 0.5d0*(ax0 + ax1)*dt
vy = vy + 0.5d0*(ay0 + ay1)*dt
vz = vz + 0.5d0*(az0 + az1)*dt

!
! update global time
!
global_time = global_time + dt
t(:) = t(:) + dt

!
! get timesteps
!
dt_old = dt
do i = 1, N
   a1 = sqrt( ax1(i)*ax1(i) + ay1(i)*ay1(i) + az1(i)*az1(i) )
   a1dot = sqrt( axdot1(i)*axdot1(i) + aydot1(i)*aydot1(i) + azdot1(i)*azdot1(i) )
   dt_try(i) = eta_dt * a1/a1dot
end do
dt = minval(dt_try)
dt = min(dt, dt_max)
if (dt .lt. dt_min) then
   print *, "VERLET: preventing dt from falling below dt_min", dt, dt_min
   dt = dt_min
end if
! prevent going past output time...
if ( (global_time + dt) .gt. time_at_next_snapshot ) then
   dt = time_at_next_snapshot - global_time
   print *, "VERLET: dt limited by snapshot time", global_time, dt
end if

end subroutine verlet_2ndorder_global_dt 
