program pos_vel_generation

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  double precision, parameter :: pi = 3.14159265358979323846_dp

  integer, parameter :: N = 70
  integer, parameter :: N_div = 2
  double precision, parameter :: D = 2.6

  real, parameter :: sigma = 1.0

  integer :: i, j, k

  integer :: nstars
  integer :: gen
  integer :: max_stars

  real :: prob_parent
  real :: dx, dy, dz

  real :: dist
  double precision :: rad_cluster
  real :: sigma_gen

  real :: random_noise

  real, dimension(:), allocatable :: x, y, z
  real, dimension(:), allocatable :: vx, vy, vz

  logical, dimension(:), allocatable :: is_parent


!
!allocate arrays
!
  allocate(x(N))
  allocate(y(N))
  allocate(z(N))

  allocate(vx(N))
  allocate(vy(N))
  allocate(vz(N))

  allocate(is_parent(N))


  max_stars = N * N_div**3

!
!first parent
!
  x(1) = 0.5
  y(1) = 0.5
  z(1) = 0.5

  vx(1) = gaussian_random(sigma)
  vy(1) = gaussian_random(sigma)
  vz(1) = gaussian_random(sigma)

  is_parent(1) = .true.
  nstars = 1
  gen = 1


  do while (nstars < max_stars)

    call generate_child(x, y, z, vx, vy, vz, is_parent, nstars, gen, N, N_div, sigma, D)

    gen  = gen + 1

  end do

  rad_cluster = 0.5
  call pruning_sphere(x, y, z, vx, vy, vz, nstars, N, rad_cluster)

  call pruning_excess(x, y, z, vx, vy, vz, nstars, N)


!
! plotting
!
  open(1, file = 'positions_velocity.dat', status = 'replace')
  
    do i = 1, nstars
      
      write(1,*) x(i), y(i), z(i), vx(i), vy(i), vz(i)

    end do
  close(1)



  deallocate(x, y, z, vx, vy, vz, is_parent)





contains
!
! gaussian random value generator
!
  real function gaussian_random(sigma)

    real, intent(in) :: sigma
    real :: u1, u2

    call random_number(u1)
    call random_number(u2)
    gaussian_random = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)

  end function gaussian_random

!
! random number generator
!
  real function rand()
    real :: rand_val
    call random_number(rand_val)
    rand = rand_val
  end function rand





end program pos_vel_generation