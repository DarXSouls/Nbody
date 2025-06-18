subroutine get_accel_for_active(N, N_act, i_active, x, y, z, vx, vy, vz, m, soft, ax, ay, az, axdot, aydot, azdot, radii)
!
!
! This subroutine returns the accelations, based on the current position,
! velocity, mass, and softening.
!
! NOTE: This is the version for the block timesteps, returning updated
! accels only for the particles listed. The other accels are left as they are. 
!
! Note that the variables here are *local* (hence the argument call), as this
! subroutine needs to be called several times within the loop.
!
! This is the 'faster' version: the if statement check for i = j has been removed
! and been replaced with a smarter i loop; we've also removed much of the duplicated
! effort in the loop.
!

use global_parameters
use timing_info
use accel_loop_list

use omp_lib
use openacc

implicit none

!
! declarations of counters, etc.
!
integer :: i, j, nn, mm

!
!
! temporary variables
double precision :: dx, dy, dz, dvx, dvy, dvz
double precision :: vdotr, rad_equiv, recip_rad_equiv
double precision :: acc_denom, acc_dot_denom
double precision :: acc_const, acc_dot_const

!
! input data declarations
!
integer :: N, N_act
integer :: i_active(N) ! this is the list of active particles (those being moved)
double precision :: x(1:N), y(1:N), z(1:N)
double precision :: vx(1:N), vy(1:N), vz(1:N)
double precision :: m(N)
double precision :: soft ! only global softening for just now... need to change this

!
! output data declarations
!
double precision :: ax(1:N), ay(1:N), az(1:N)
double precision :: axdot(1:N), aydot(1:N), azdot(1:N)
double precision :: radii(1:N, 1:N)

!
! count total calls to subroutine
!
calls_to_get_accel = calls_to_get_accel + 1

!
! sanity check!
!
if (N_act == 0) then
   print *, 'No active particles to advance!'
   stop
end if

!
! loop around the active particles to get the accels etc.
!

!$omp parallel do
do nn = 1, N_act
   i = i_active(nn)
   !
   ! clear the accel arrays for the *active* particles. The other accels are left as they
   ! are (these particles are being predicted and so need to retain their old accel values)
   !   
   ax(i) = 0
   ay(i) = 0
   az(i) = 0
   axdot(i) = 0
   aydot(i) = 0
   azdot(i) = 0
   !
   ! loop around other particles
   !
   do mm = 1, N-1
     j = indices(i+mm)
     dx = x(j) - x(i)
     dy = y(j) - y(i)
     dz = z(j) - z(i)
     dvx = vx(j) - vx(i)
     dvy = vy(j) - vy(i)
     dvz = vz(j) - vz(i)
     vdotr = dx*dvx + dy*dvy + dz*dvz
     rad_equiv = sqrt(dx*dx + dy*dy + dz*dz + soft*soft)
     recip_rad_equiv = 1.0D0 / rad_equiv
     acc_denom = recip_rad_equiv * recip_rad_equiv * recip_rad_equiv
     acc_dot_denom = acc_denom * recip_rad_equiv * recip_rad_equiv
     acc_const = G_int*m(j)*acc_denom
     acc_dot_const = 3.0D0*vdotr*G_int*m(j)*acc_dot_denom

     !
     ! accel of body j on i
     !
     ax(i) = ax(i) + acc_const * dx 
     ay(i) = ay(i) + acc_const * dy
     az(i) = az(i) + acc_const * dz

     !
     ! time derivative of accel of body j on i
     !
     axdot(i) = axdot(i) + acc_const*dvx + acc_dot_const*dx
     aydot(i) = aydot(i) + acc_const*dvy + acc_dot_const*dy
     azdot(i) = azdot(i) + acc_const*dvz + acc_dot_const*dz
     !
     ! populate radii storage array
     !
     radii(i, j) = rad_equiv
     !
     ! store number of times through loop for diagnostics
     !
     accel_function_counter = accel_function_counter + 1 
  end do
end do
!$omp end parallel do

end subroutine get_accel_for_active
