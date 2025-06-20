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

use global_parameters
use timing_info

implicit none

!
! declarations of counters, etc.
!
integer :: i, j, nn

!
!
! temporary variables
double precision :: dx, dy, dz, dvx, dvy, dvz
double precision :: vdotr, rad_equiv
double precision :: acc_denom, acc_dot_denom

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
   do j = 1,  N
     if ( i .eq. j ) cycle
     dx = x(j) - x(i)
     dy = y(j) - y(i)
     dz = z(j) - z(i)
     dvx = vx(j) - vx(i)
     dvy = vy(j) - vy(i)
     dvz = vz(j) - vz(i)
     vdotr = dx*dvx + dy*dvy + dz*dvz
     rad_equiv = sqrt(dx*dx + dy*dy + dz*dz + soft*soft)
     acc_denom = 1.0D0 / rad_equiv**3
     acc_dot_denom = 1.0D0 / rad_equiv**5

     !
     ! accel of body j on i
     !
     ax(i) = ax(i) + G_int*m(j)*dx*acc_denom
     ay(i) = ay(i) + G_int*m(j)*dy*acc_denom
     az(i) = az(i) + G_int*m(j)*dz*acc_denom

     !
     ! time derivative of accel of body j on i
     !
     axdot(i) = axdot(i) + G_int*m(j)*(dvx*acc_denom + 3.0*vdotr*dx*acc_dot_denom)
     aydot(i) = aydot(i) + G_int*m(j)*(dvy*acc_denom + 3.0*vdotr*dy*acc_dot_denom)
     azdot(i) = azdot(i) + G_int*m(j)*(dvz*acc_denom + 3.0*vdotr*dz*acc_dot_denom)
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
