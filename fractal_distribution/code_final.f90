program particle_generation
    implicit none

    integer, parameter :: N = 70, Ndiv = 2
    real, parameter :: D = 2.6, sigma = 1.0

    integer :: i, j, k, nstars, generation
    real :: prob_parent, dx, dy, dz, dist, sigma_g
    
    real, dimension(:), allocatable :: x, y, z, vx, vy, vz, m
    real, dimension(:), allocatable :: mass_samples
    logical, dimension(:), allocatable :: is_parent

    integer :: max_stars
    integer :: sample_index

    allocate(x(N*8), y(N*8), z(N*8))
    allocate(vx(N*8), vy(N*8), vz(N*8))
    allocate(m(N*8))
    allocate(mass_samples(1000))
    allocate(is_parent(N*8))

    open(10, file='mass_sample1.dat', status='old')
    do i = 1, 1000
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
        call generate_children(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        generation = generation + 1
    end do

    ! Remove parents
    call remove_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)

    ! Prune to sphere
    call prune_to_sphere(x, y, z, vx, vy, vz, m, nstars)

    ! Prune excess stars
    call remove_excess_stars(x, y, z, vx, vy, vz, m, nstars, N)

    ! Scaling velocities for virial ratio
    call scale_velocities(vx, vy, vz, m, nstars, 0.5)  ! Assuming alpha = 0.5

    ! Output positions and velocities
    open(1, file='positions_velocity.dat', status='replace')
    do i = 1, N
        write(1,*) x(i), y(i), z(i), vx(i), vy(i), vz(i), m(i)
    end do
    close(1)

    deallocate(x, y, z, vx, vy, vz, m, is_parent, mass_samples)

contains

    subroutine generate_children(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        integer, intent(inout) :: nstars, generation, sample_index
        integer, intent(in) :: N, Ndiv
        real, intent(in) :: sigma, D
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        logical, dimension(:), intent(inout) :: is_parent
        integer :: i, j, k, parent_idx, child_idx, new_nstars
        real :: dx, dy, dz, prob_parent, sigma_g

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

                            dx = real(i) / Ndiv
                            dy = real(j) / Ndiv
                            dz = real(k) / Ndiv
                            x(child_idx) = x(parent_idx) + dx
                            y(child_idx) = y(parent_idx) + dy
                            z(child_idx) = z(parent_idx) + dz 

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
    end subroutine generate_children

    subroutine remove_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)
        integer, intent(inout) :: nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        logical, dimension(:), intent(inout) :: is_parent
        integer :: i, j
        integer, dimension(:), allocatable :: to_remove
        integer :: remove_count

        allocate(to_remove(nstars))
        remove_count = 0

        do i = 1, nstars
            if (is_parent(i)) then
                remove_count = remove_count + 1
                to_remove(remove_count) = i
            end if
        end do

        do i = remove_count, 1, -1
            do j = to_remove(i), nstars-1
                x(j) = x(j+1)
                y(j) = y(j+1)
                z(j) = z(j+1)
                vx(j) = vx(j+1)
                vy(j) = vy(j+1)
                vz(j) = vz(j+1)
                m(j) = m(j+1)
            end do
            nstars = nstars - 1
        end do

        deallocate(to_remove)
    end subroutine remove_parents

    subroutine prune_to_sphere(x, y, z, vx, vy, vz, m, nstars)
        integer, intent(inout) :: nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        integer :: i, j
        real :: dist
        integer, dimension(:), allocatable :: to_remove
        integer :: remove_count

        allocate(to_remove(nstars))
        remove_count = 0

        do i = 1, nstars
            dist = sqrt((x(i)-0.5)**2 + (y(i)-0.5)**2 + (z(i)-0.5)**2)
            if (dist > 0.5) then  ! 0.5 is the radius of the sphere
                remove_count = remove_count + 1
                to_remove(remove_count) = i
            end if
        end do

        do i = remove_count, 1, -1
            do j = to_remove(i), nstars-1
                x(j) = x(j+1)
                y(j) = y(j+1)
                z(j) = z(j+1)
                vx(j) = vx(j+1)
                vy(j) = vy(j+1)
                vz(j) = vz(j+1)
                m(j) = m(j+1)
            end do
            nstars = nstars - 1
        end do

        deallocate(to_remove)
    end subroutine prune_to_sphere

    subroutine remove_excess_stars(x, y, z, vx, vy, vz, m, nstars, target_nstars)
        integer, intent(inout) :: nstars
        integer, intent(in) :: target_nstars
        real, dimension(:), intent(inout) :: x, y, z, vx, vy, vz, m
        integer :: i, j, idx
        real :: rnd

        do i = 1, nstars-target_nstars
            call random_number(rnd)
            idx = int(rnd * nstars) + 1
            do j = idx, nstars-1
                x(j) = x(j+1)
                y(j) = y(j+1)
                z(j) = z(j+1)
                vx(j) = vx(j+1)
                vy(j) = vy(j+1)
                vz(j) = vz(j+1)
                m(j) = m(j+1)
            end do
            nstars = nstars - 1
        end do
    end subroutine remove_excess_stars

    subroutine scale_velocities(vx, vy, vz, m, nstars, alpha)
        real, intent(inout) :: vx(:), vy(:), vz(:), m(:)
        integer, intent(in) :: nstars
        real, intent(in) :: alpha
        real :: T, V, scale_factor
        integer :: i, j

        T = 0.0
        V = 0.0

        ! Calculate total kinetic energy T
        do i = 1, nstars
            T = T + 0.5 * m(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2)
        end do

        ! Calculate total potential energy V (simple approximation assuming uniform mass distribution)
        do i = 1, nstars
            do j = i+1, nstars
                V = V - m(i) * m(j) / sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
            end do
        end do

        scale_factor = sqrt((alpha * abs(V)) / T)

        ! Scale velocities
        do i = 1, nstars
            vx(i) = vx(i) * scale_factor
            vy(i) = vy(i) * scale_factor
            vz(i) = vz(i) * scale_factor
        end do
    end subroutine scale_velocities

    real function gaussian_random(sigma)
        real, intent(in) :: sigma
        real :: u1, u2

        call random_number(u1)
        call random_number(u2)
        gaussian_random = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * 3.141592653589793 * u2)
    end function gaussian_random

    real function rand()
        real :: rand_val
        call random_number(rand_val)
        rand = rand_val
    end function rand

end program particle_generation
