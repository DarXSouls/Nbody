program particle_generation
    implicit none

    ! Parameters
    integer, parameter :: N = 70
    integer, parameter :: Ndiv = 2
    real, parameter :: D = 2.6
    real :: vir = 0.5

    ! Standard deviation for Gaussian velocity
    real, parameter :: sigma = 1.0

    integer :: i, j, k
    integer :: nstars, generation
    integer :: sample_index
    real :: prob_parent, sigma_g
    real :: dx, dy, dz, dist
    real :: vir_frac
    real :: vcom_x, vcom_y, vcom_z
    real :: total_mass

    real, dimension(:), allocatable :: x, y, z
    real, dimension(:), allocatable :: vx, vy, vz
    real, dimension(:), allocatable :: m
    real, dimension(:), allocatable :: mass_samples
    logical, dimension(:), allocatable :: is_parent

    integer :: max_stars

    ! Physical constants
    double precision, parameter :: gg = 6.672d-08
    double precision, parameter :: year = 3600.*24.*365.2425
    double precision, parameter :: solar_mass = 1.98847d33
    double precision, parameter :: earth_mass = 5.9722d27
    double precision, parameter :: au = 1.495978707d13
    double precision, parameter :: pc = 3.085677581d18

    ! Units
    double precision :: units_mass_g
    double precision :: units_length_cm
    double precision :: units_time_s
    real :: G_int

    units_mass_g = solar_mass
    units_length_cm = pc
    units_time_s = year
    G_int = gg * units_mass_g * units_time_s**2 / units_length_cm**3

    max_stars = N * Ndiv**3
    
    allocate(x(N*8), y(N*8), z(N*8))
    allocate(vx(N*8), vy(N*8), vz(N*8))
    allocate(m(N*8))
    allocate(mass_samples(1000000))
    allocate(is_parent(N*8))

    ! Open the mass sample file
    open(10, file='mass_sample2.dat', status='old', iostat=i)
    if (i /= 0) then
        print *, "Error: Cannot open file 'mass_sample2.dat'"
        stop
    end if
    do i = 1, 1000000
        read(10,*) mass_samples(i)
    end do
    close(10)

    sample_index = 1
    max_stars = N * Ndiv**3

    x(1) = 0.5
    y(1) = 0.5
    z(1) = 0.5
    vx(1) = gaussian_random(sigma)
    vy(1) = gaussian_random(sigma)
    vz(1) = gaussian_random(sigma)

    m(1) = mass_samples(sample_index)
    sample_index = sample_index + 1

    is_parent(1) = .true.
    nstars = 1
    generation = 1

    ! Generate children
    do while (nstars < N)
        call gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        generation = generation + 1
    end do

    ! Remove parents
    call prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)

    ! Prune sphere
    call prune_sphere(x, y, z, vx, vy, vz, m, nstars)

    ! Prune excess
    call prune_excess(x, y, z, vx, vy, vz, m, nstars, N)

    ! Subtract the center of mass velocity
    call subtract_com_velocity(vx, vy, vz, m, nstars)

    ! Scale velocities
    call scale_vel(vx, vy, vz, m, nstars, vir, G_int)

    ! File writing: pos, vel, mass
    open(1, file='pos_vel_m2.dat', status='replace')
    do i = 1, N
        write(1, 4000) x(i), y(i), z(i), vx(i), vy(i), vz(i), m(i)
    end do
    4000 format(16(E20.14, 1X))
    close(1)

    deallocate(x, y, z, vx, vy, vz, m, is_parent, mass_samples)

contains

    ! Subroutine to generate child particles
    subroutine gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        implicit none
        integer, intent(inout) :: nstars, generation, sample_index
        integer, intent(in) :: N, Ndiv
        real, intent(in) :: D, sigma
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        logical, dimension(:), intent(inout) :: is_parent

        integer :: i, j, k, parent_idx, child_idx, new_nstars
        real :: dx, dy, dz, prob_parent, sigma_g
        real :: noise = 0.05

        logical :: is_duplicate

        new_nstars = nstars
        sigma_g = sigma * (1.0 / Ndiv) ** generation

        do parent_idx = 1, nstars
            if (is_parent(parent_idx)) then
                do i = 0, Ndiv-1
                    do j = 0, Ndiv-1
                        do k = 0, Ndiv-1
                            if (new_nstars >= N) exit
                            new_nstars = new_nstars + 1
                            child_idx = new_nstars

                            dx = (i + 0.5) / Ndiv
                            dy = (j + 0.5) / Ndiv
                            dz = (k + 0.5) / Ndiv

                            x(child_idx) = x(parent_idx) + dx + gaussian_random(noise)
                            y(child_idx) = y(parent_idx) + dy + gaussian_random(noise)
                            z(child_idx) = z(parent_idx) + dz + gaussian_random(noise)

                            vx(child_idx) = vx(parent_idx) + gaussian_random(sigma_g)
                            vy(child_idx) = vy(parent_idx) + gaussian_random(sigma_g)
                            vz(child_idx) = vz(parent_idx) + gaussian_random(sigma_g)

                            m(child_idx) = mass_samples(sample_index)
                            sample_index = sample_index + 1

                            prob_parent = real(Ndiv) ** (D - 3)
                            is_parent(child_idx) = (rand() < prob_parent)
                            

                        end do
                    end do
                end do
                is_parent(parent_idx) = .false.
            end if
        end do

        nstars = new_nstars
    end subroutine gen_child

    ! Subroutine to subtract the center of mass velocity
    subroutine subtract_com_velocity(vx, vy, vz, m, nstars)
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(in) :: m
        integer, intent(in) :: nstars

        integer :: i
        real :: vcom_x, vcom_y, vcom_z
        real :: total_mass

        vcom_x = 0.0
        vcom_y = 0.0
        vcom_z = 0.0
        total_mass = 0.0

        ! Calculate center of mass velocity
        do i = 1, nstars
            vcom_x = vcom_x + m(i) * vx(i)
            vcom_y = vcom_y + m(i) * vy(i)
            vcom_z = vcom_z + m(i) * vz(i)
            total_mass = total_mass + m(i)
        end do

        vcom_x = vcom_x / total_mass
        vcom_y = vcom_y / total_mass
        vcom_z = vcom_z / total_mass

        ! Subtract it from each particle's velocity
        do i = 1, nstars
            vx(i) = vx(i) - vcom_x
            vy(i) = vy(i) - vcom_y
            vz(i) = vz(i) - vcom_z
        end do
    end subroutine subtract_com_velocity

    ! Subroutine to prune parent particles
    subroutine prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)
        integer, intent(inout) :: nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        logical, dimension(:), intent(inout) :: is_parent

        integer :: i, j

        do i = nstars, 1, -1
            if (is_parent(i)) then
                do j = i, nstars - 1
                    x(j) = x(j+1)
                    y(j) = y(j+1)
                    z(j) = z(j+1)
                    vx(j) = vx(j+1)
                    vy(j) = vy(j+1)
                    vz(j) = vz(j+1)
                    m(j) = m(j+1)
                end do
                nstars = nstars - 1
            end if
        end do
    end subroutine prune_parents

    ! Subroutine to prune particles outside a sphere
    subroutine prune_sphere(x, y, z, vx, vy, vz, m, nstars)
        integer, intent(inout) :: nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m

        integer :: i, j
        real :: dist, radius

        radius = 1.0
        do i = nstars, 1, -1
            dist = sqrt(x(i)**2 + y(i)**2 + z(i)**2)
            if (dist > radius) then
                do j = i, nstars - 1
                    x(j) = x(j+1)
                    y(j) = y(j+1)
                    z(j) = z(j+1)
                    vx(j) = vx(j+1)
                    vy(j) = vy(j+1)
                    vz(j) = vz(j+1)
                    m(j) = m(j+1)
                end do
                nstars = nstars - 1
            end if
        end do
    end subroutine prune_sphere

    ! Subroutine to prune excess particles
    subroutine prune_excess(x, y, z, vx, vy, vz, m, nstars, N)
        integer, intent(inout) :: nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        integer, intent(in) :: N

        if (nstars > N) then
            nstars = N
        end if
    end subroutine prune_excess

    ! Subroutine to scale velocities
    subroutine scale_vel(vx, vy, vz, m, nstars, vir, G_int)
        integer, intent(in) :: nstars
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(in) :: m
        real, intent(in) :: vir, G_int

        integer :: i, j
        real :: ke, pe, scaling_factor
        real :: dist

        ke = 0.0
        pe = 0.0

        ! Calculate kinetic energy
        do i = 1, nstars
            ke = ke + 0.5 * m(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2)
        end do

        ! Calculate potential energy
        do i = 1, nstars
            do j = i+1, nstars
                dist = sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
                pe = pe - G_int * m(i) * m(j) / dist
            end do
        end do

        ! Scale velocities
        scaling_factor = sqrt(2.0 * vir * abs(pe) / ke)
        do i = 1, nstars
            vx(i) = vx(i) * scaling_factor
            vy(i) = vy(i) * scaling_factor
            vz(i) = vz(i) * scaling_factor
        end do
    end subroutine scale_vel

    ! Function to generate Gaussian random numbers
    real function gaussian_random(sigma)
        real, intent(in) :: sigma
        real :: u1, u2, z0

        u1 = rand()
        u2 = rand()

        z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
        gaussian_random = z0 * sigma
    end function gaussian_random
end program particle_generation
