subroutine generate_child(x, y, z, vx, vy, vz, is_parent, nstars, gen, N, N_div, sigma, D, max_stars)

implicit none


integer, parameter :: dp = kind(1.0d0)
double precision, parameter :: pi = 3.14159265358979323846_dp

integer, intent(inout) :: nstars, gen

integer, intent(in) :: N 
integer, intent(in) :: N_div 
double precision, intent(in) :: D 
integer, intent(in) :: max_stars

real, intent(in) :: sigma

real, dimension(:), intent(inout) :: x, y, z
real, dimension(:), intent(inout) :: vx, vy, vz

logical, dimension(:), intent(inout) :: is_parent

integer :: i, j, k, l
integer :: parent_idx, child_idx
integer :: new_nstars

real :: dx, dy, dz
real :: prob_parent
real :: sigma_gen



new_nstars = nstars
sigma_gen = sigma * (1.0 / N_div) ** gen

do parent_idx = 1, nstars
   if (is_parent(parent_idx)) then
       do i = 0, N_div - 1
           do j = 0, N_div - 1
               do k = 0, N_div - 1

                  
                   if (new_nstars >= max_stars) exit

                    new_nstars = new_nstars + 1
                    child_idx = new_nstars

                    dx = real(i) / N_div
                    dy = real(j) / N_div
                    dz = real(k) / N_div

                    !pos  
                    x(child_idx) = x(parent_idx) + dx
                    y(child_idx) = y(parent_idx) + dy
                    z(child_idx) = z(parent_idx) + dz            
                    
                    !vel
                    vx(child_idx) = vx(parent_idx) + gaussian_random(sigma_gen)
                    vy(child_idx) = vy(parent_idx) + gaussian_random(sigma_gen)
                    vz(child_idx) = vz(parent_idx) + gaussian_random(sigma_gen)
                    
                    prob_parent = (N_div ** (D - 3))
                    is_parent(child_idx) = (rand() < prob_parent)

                end do
            end do
        end do
        is_parent(parent_idx) = .false.
    end if
end do
    
        
nstars = new_nstars



contains
!
!gaussian random value generator
!
real function gaussian_random(sigma)

real, intent(in) :: sigma
real :: u1, u2

call random_number(u1)
call random_number(u2)
gaussian_random = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)

end function gaussian_random



!
!random number generator
!
real function rand()
real :: rand_val
call random_number(rand_val)
rand = rand_val
end function rand



end subroutine generate_child