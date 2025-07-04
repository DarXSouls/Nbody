subroutine get_accel_info(N, x, y, z, vx, vy, vz, m, soft, ax, ay, az, axdot, aydot, azdot)
!
!
! This subroutine returns the accelations, based on the current position,
! velocity, mass, and softening.
!
! NOTE: This is the version for *global timesteps*, and has the slow, full N^2
! loop (i.e. does not make use of the symmetry in the force equation). It gets
! the accels and their derivatives for all particles at once.
!
! Note that the variables here are *local* (hence the argument call), as this
! subroutine needs to be called several times within the loop.
!

use global_parameters

use omp_lib
use openacc

implicit none

!
! declarations of counters, etc.
!
integer :: i, j

!
!
! temporary variables
double precision :: dx, dy, dz, dvx, dvy, dvz
double precision :: vdotr, rad_equiv
double precision :: acc_denom, acc_dot_denom

!
! input data declarations
!
integer N
double precision :: x(1:N), y(1:N), z(1:N)
double precision :: vx(1:N), vy(1:N), vz(1:N)
double precision :: m(N)
double precision :: soft ! only global softening for just now... need to change this

!
! output data declarations
!
double precision :: ax(1:N), ay(1:N), az(1:N)
double precision :: axdot(1:N), aydot(1:N), azdot(1:N)

!
! clear the accel arrays
!
ax(:) = 0
ay(:) = 0
az(:) = 0
axdot(:) = 0
aydot(:) = 0
azdot(:) = 0

!
! N^2 loop around the particles to get the accels etc.
!

!$omp parallel do
do i = 1, N
   do j = 1, N
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
     axdot(i) = axdot(i) + G_int*m(j)*(dvx*acc_denom + 3.0d0*vdotr*dx*acc_dot_denom)
     aydot(i) = aydot(i) + G_int*m(j)*(dvy*acc_denom + 3.0d0*vdotr*dy*acc_dot_denom)
     azdot(i) = azdot(i) + G_int*m(j)*(dvz*acc_denom + 3.0d0*vdotr*dz*acc_dot_denom)
  end do
end do
!$omp end parallel do


end subroutine get_accel_info
