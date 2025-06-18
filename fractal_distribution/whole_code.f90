subroutine particle_generation
    
    use global_parameters
    use global_variables
    
    implicit none

    ! params
    
    integer, parameter :: Ndiv = 2.0
    real(8), parameter :: D = 2.6
    real :: vir = 2.0

    ! std for gauss vel
    real(8), parameter :: sigma = 1.0

    integer :: i, j, k
    integer :: nstars, generation
    real :: prob_parent, sigma_g
    real :: dx, dy, dz, dist
    
    !real, dimension(:), allocatable :: x, y, z
    !real, dimension(:), allocatable :: vx, vy, vz
    !real, dimension(:), allocatable :: m
    real, dimension(:), allocatable :: mass_samples
    logical, dimension(:), allocatable :: is_parent

    integer :: max_stars
    integer :: sample_index

    real :: vir_frac



    real(8) :: Time_cross
    real(8) :: Time_evol
    real(8) :: rad_cluster = 0.5

    integer :: rand_star
    real :: rnd_no


    !
    ! Physical constants
    !
    !double precision, parameter :: gg = 6.672d-08
    !double precision, parameter :: year = 3600.*24.*365.2425
    !double precision, parameter :: solar_mass = 1.98847d33
    double precision, parameter :: earth_mass = 5.9722d27
    double precision, parameter :: au = 1.495978707d13
    double precision, parameter :: pc = 3.085677581d18
    !
    ! set our units
    !
    !double precision :: units_mass_g
    !double precision :: units_length_cm
    !double precision :: units_time_s
    !double precision :: G_int

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
    open(10, file='mass_sample6.dat', status='old')
    do i = 1, 10000
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
    !vx(1) = 0.0
    !vy(1) = 0.0
    !vz(1) = 0.0

    m(1) = mass_samples(sample_index)
    sample_index = sample_index + 1
    
    is_parent(1) = .true.
    nstars = 1
    generation = 1


    ! gen children
    do while (nstars < N*4)
        call gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)
        generation = generation + 1
    end do

    ! remove parent
    call prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)

    ! Prune sphere
    call prune_sphere(x, y, z, vx, vy, vz, m, nstars, rad_cluster)

    ! prune excess
    call prune_excess(x, y, z, vx, vy, vz, m, nstars, N)

    ! vcom substraction
    call vcom_subs(vx, vy, vz, m, nstars)

    ! scaling vel
    call scale_vel(vx, vy, vz,  m, nstars, vir, real(G_int, kind=4))


    !call random_number(rnd_no)
    !rand_star = int(rnd_no * nstars) + 1
    !m(rand_star) = 10.0


    

    !
    ! Crossing times
    !
    print *, 'total mass: ', sum(m(1:N))

    Time_cross = (sqrt(vir / (4*G_int*sum(m(1:N)))) * (rad_cluster)**(3/2))
    print *, 'Crossing time: ', Time_cross

    Time_evol = Time_cross * (0.1*N / log(real(N))) * 100
    print *, 'Evolution time: ', Time_evol


    !
    ! file writing pos vel m
    !
    open(1, file = 'ics.dat', status = 'replace')

      !write(1, *) 'Crossing time :', Time_cross

      do i = 1, N

          write(1,4000) x(i), y(i), z(i), vx(i), vy(i), vz(i), m(i)

        end do
        4000 format(16(E20.14, 1X))

    close(1)


    
    deallocate(x, y, z, vx, vy, vz, m, is_parent, mass_samples)



contains

    subroutine gen_child(x, y, z, vx, vy, vz, m, is_parent, nstars, generation, N, Ndiv, sigma, D, sample_index)

        use global_parameters
        !use global_variables

        implicit none

        integer, intent(inout) :: nstars, generation, sample_index
        
        integer, intent(in) :: N
        integer, intent(in) :: Ndiv
        real(8), intent(in) :: D

        real(8), intent(in) :: sigma

        real(8), dimension(:), intent(inout) :: x, y, z
        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(inout) :: m

        logical, dimension(:), intent(inout) :: is_parent

        integer :: i, j, k
        integer :: parent_idx, child_idx
        integer :: new_nstars

        real(8) :: dx, dy, dz
        real(8) :: prob_parent, sigma_g

        ! noise
        real(8) :: noise
        noise = 0.05

        new_nstars = nstars
        sigma_g = sigma * (1.0 / Ndiv) ** generation

        do parent_idx = 1, nstars
            if (is_parent(parent_idx)) then
                do i = 0, Ndiv-1
                    do j = 0, Ndiv-1
                        do k = 0, Ndiv-1
                            !if (new_nstars >= N) exit
                            
            
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
                            
                            !end if

                            
                        end do
                    end do
                end do
                is_parent(parent_idx) = .false.
            end if
        end do

        nstars = new_nstars
    end subroutine gen_child



    subroutine prune_parents(x, y, z, vx, vy, vz, m, is_parent, nstars)

        use global_parameters
        !use global_variables

        implicit none
        
        integer, intent(inout) :: nstars

        real(8), dimension(:), intent(inout) :: x, y, z
        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(inout) :: m

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
    subroutine prune_sphere(x, y, z, vx, vy, vz, m, nstars, rad_cluster)

        use global_parameters
        !use global_variables

        implicit none
        
        integer, intent(inout) :: nstars

        real(8), dimension(:), intent(inout) :: x, y, z
        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(inout) :: m

        integer :: i, j
        real :: dist
        real(8), intent(inout) :: rad_cluster

        integer, dimension(:), allocatable :: to_remove
        integer :: remove_count

        allocate(to_remove(nstars))
        remove_count = 0

        

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

        use global_parameters
        !use global_variables

        implicit none
        
        integer, intent(inout) :: nstars
        integer, intent(in) :: target_nstars

        real(8), dimension(:), intent(inout) :: x, y, z
        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(inout) :: m

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

        use global_parameters
        !use global_variables

        implicit none
        
        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(in) :: m
        integer, intent(in) :: nstars

        integer :: i
        real :: vcom_x, vcom_y, vcom_z
        real :: total_mass

        vxcom = 0.0
        vycom = 0.0
        vzcom = 0.0
        total_mass = 0.0

        do i = 1, nstars
            vxcom = vxcom + m(i) * vx(i)
            vycom = vycom + m(i) * vy(i)
            vzcom = vzcom + m(i) * vz(i)
            total_mass = total_mass + m(i)
        end do

        vxcom = vxcom / total_mass
        vycom = vycom / total_mass
        vzcom = vzcom / total_mass

        do i = 1, nstars
            vx(i) = vx(i) - vxcom
            vy(i) = vy(i) - vycom
            vz(i) = vz(i) - vzcom
        end do

    end subroutine vcom_subs



    ! scaling vels using energies
    subroutine scale_vel(vx, vy, vz, m, nstars, vir, G_int)

        !use global_parameters
        !use global_variables

        implicit none

        real(8), dimension(:), intent(inout) :: vx, vy, vz
        real(8), dimension(:), intent(inout) :: m
        integer, intent(in) :: nstars
        real, intent(in) :: vir
        real, intent(in) :: G_int

        integer :: i, j
        real :: e_kin, e_grav
        real :: scale_factor

        
        

        ! E kin
        e_kin = 0.5d0 * sum(m*(vx*vx + vy*vy + vz*vz))

        e_grav = 0.0

        ! E grav
        do i = 1, nstars - 1
            do j = i + 1, nstars
                
                e_grav = e_grav + G_int*m(i)*m(j) / sqrt( (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2 )
            end do
        end do

        print *, 'E_kin: ', e_kin
        print *, 'E_grav: ', e_grav
        print *, 'e_rat :', vir, e_grav/e_kin


        scale_factor = sqrt((abs(e_grav)/e_kin)/vir)
        print *, 'vir frac: ', scale_factor

        ! scaling yeahboi
        vx = vx * scale_factor
        vy = vy * scale_factor
        vz = vz * scale_factor

        e_kin = 0.5d0 * sum(m*(vx*vx + vy*vy + vz*vz))

        print *, 'new e_kin: ', e_kin

    end subroutine scale_vel



    real function gaussian_random(sigma)
        implicit none 

        real(8), intent(in) :: sigma
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




end subroutine particle_generation
