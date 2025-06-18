program particle_generation
    
    use global_parameters
    use global_variables
    
    implicit none

    ! params
    
    integer, parameter :: Ndiv = 2
    real, parameter :: D = 2.6
    real :: vir = 2.5

    ! std for gauss vel
    real, parameter :: sigma = 1.0

    integer :: i, j, k
    integer :: nstars, generation
    real :: prob_parent, sigma_g
    real :: dx, dy, dz, dist
    
    real, dimension(:), allocatable :: x, y, z
    real, dimension(:), allocatable :: vx, vy, vz
    real, dimension(:), allocatable :: m
    real, dimension(:), allocatable :: mass_samples
    logical, dimension(:), allocatable :: is_parent

    integer :: max_stars
    integer :: sample_index

    real :: vir_frac

    !
    ! Physical constants
    !
    double precision, parameter :: gg = 6.672d-08
    double precision, parameter :: year = 3600.*24.*365.2425
    double precision, parameter :: solar_mass = 1.98847d33
    double precision, parameter :: earth_mass = 5.9722d27
    double precision, parameter :: au = 1.495978707d13
    double precision, parameter :: pc = 3.085677581d18
    !
    ! set our units
    !
    double precision :: units_mass_g
    double precision :: units_length_cm
    double precision :: units_time_s
    double precision :: G_int

    units_mass_g = solar_mass
    units_length_cm = pc
    units_time_s = year
    G_int = gg * units_mass_g * units_time_s**2 / units_length_cm**3

    
    allocate(x(N*8), y(N*8), z(N*8))
    allocate(vx(N*8), vy(N*8), vz(N*8))

    allocate(m(N*8))
    allocate(mass_samples(1000000))
    
    allocate(is_parent(N*8))

    ! opening the mass sample file
    open(10, file='mass_sample3.dat', status='old')
    do i = 1, 500
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


    ! gen children
    do while (nstars < N*8)
        call gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        generation = generation + 1
    end do

    ! remove parent
    call prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)

    ! Prune sphere
    call prune_sphere(x, y, z, vx, vy, vz, m, nstars)

    ! prune excess
    call prune_excess(x, y, z, vx, vy, vz, m, nstars, N)

    ! vcom substraction
    call vcom_subs(vx, vy, vz, m, nstars)

    ! scaling vel
    call scale_vel(vx, vy, vz,  m, nstars, vir, real(G_int, kind=4))


    !
    ! file writing pos vel m
    !
    open(1, file = 'pos_vel_m.dat', status = 'replace')

      do i = 1, N

          write(1,4000) x(i), y(i), z(i), vx(i), vy(i), vz(i), m(i)

        end do
        4000 format(16(E20.14, 1X))

    close(1)


    
    deallocate(x, y, z, vx, vy, vz, m, is_parent, mass_samples)



contains

    subroutine gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        integer, intent(inout) :: nstars, generation, sample_index
        
        integer, intent(in) :: N
        integer, intent(in) :: Ndiv
        real, intent(in) :: D

        real, intent(in) :: sigma

        real, dimension(:), intent(inout) :: x, y, z
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(inout) :: m

        logical, dimension(:), intent(inout) :: is_parent

        integer :: i, j, k
        integer :: parent_idx, child_idx
        integer :: new_nstars

        real :: dx, dy, dz
        real :: prob_parent, sigma_g

        ! noise
        real :: noise
        noise = 0.05

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



    subroutine prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)
        
        integer, intent(inout) :: nstars

        real, dimension(:), intent(inout) :: x, y, z
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(inout) :: m

        logical, dimension(:), intent(inout) :: is_parent
        integer, dimension(:), allocatable :: to_remove

        integer :: i, j
        integer :: remove_count

        allocate(to_remove(nstars))
        remove_count = 0

        ! which parents pooof
        do i = 1, nstars
            if (is_parent(i)) then
                remove_count = remove_count + 1
                to_remove(remove_count) = i
            end if
        end do

        ! remaining stars shifting down
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
    end subroutine prune_parents



    ! removing the stars outside the sphere
    subroutine prune_sphere(x, y, z, vx, vy, vz, m, nstars)
        
        integer, intent(inout) :: nstars

        real, dimension(:), intent(inout) :: x, y, z
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(inout) :: m

        integer :: i, j
        real :: dist
        real :: rad_cluster

        integer, dimension(:), allocatable :: to_remove
        integer :: remove_count

        allocate(to_remove(nstars))
        remove_count = 0

        rad_cluster = 0.5

        do i = 1, nstars
            dist = sqrt((x(i)-0.5)**2 + (y(i)-0.5)**2 + (z(i)-0.5)**2)
            if (dist > rad_cluster) then  
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
    end subroutine prune_sphere



    ! removing excess
    subroutine prune_excess(x, y, z, vx, vy, vz, m, nstars, target_nstars)
        
        integer, intent(inout) :: nstars
        integer, intent(in) :: target_nstars

        real, dimension(:), intent(inout) :: x, y, z
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(inout) :: m

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
    end subroutine prune_excess



    ! vcom substraction
    subroutine vcom_subs(vx, vy, vz, m, nstars)
        
        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(in) :: m
        integer, intent(in) :: nstars

        integer :: i
        real :: vcom_x, vcom_y, vcom_z
        real :: m_total

        vcom_x = 0.0
        vcom_y = 0.0
        vcom_z = 0.0
        m_total = 0.0

        do i = 1, nstars
            vcom_x = vcom_x + m(i) * vx(i)
            vcom_y = vcom_y + m(i) * vy(i)
            vcom_z = vcom_z + m(i) * vz(i)
            m_total = m_total + m(i)
        end do

        vcom_x = vcom_x / m_total
        vcom_y = vcom_y / m_total
        vcom_z = vcom_z / m_total

        do i = 1, nstars
            vx(i) = vx(i) - vcom_x
            vy(i) = vy(i) - vcom_y
            vz(i) = vz(i) - vcom_z
        end do

    end subroutine vcom_subs



    ! scaling vels using energies
    subroutine scale_vel(vx, vy, vz, m, nstars, vir, G_int)

        real, dimension(:), intent(inout) :: vx, vy, vz
        real, dimension(:), intent(inout) :: m
        integer, intent(in) :: nstars
        real, intent(in) :: vir
        real, intent(in) :: G_int

        integer :: i, j
        real :: E_kin, E_grav
        real :: scale_factor

        E_kin = 0.0
        E_grav = 0.0

        ! E kin
        do i = 1, nstars

          E_kin = E_kin + 0.5 * m(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2)

        end do

        ! E grav
        do i = 1, nstars - 1
            do j = i + 1, nstars
                
                E_grav = E_grav - G_int * m(i) * m(j) / sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
            end do
        end do

        scale_factor = sqrt(vir * abs(E_grav) / E_kin )

        ! scaling yeahboi
        do i = 1, nstars
            
            vx(i) = vx(i) * scale_factor
            vy(i) = vy(i) * scale_factor
            vz(i) = vz(i) * scale_factor
        end do

    end subroutine scale_vel



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
