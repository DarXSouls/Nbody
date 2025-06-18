subroutine pruning_excess(x, y, z, vx, vy, vz, nstars, N)

integer, intent(inout) :: nstars
integer, intent(inout) :: N

real, dimension(:), intent(inout) :: x, y, z
real, dimension(:), intent(inout) :: vx, vy, vz

integer :: i, j
integer :: idx

do i = 1, nstars - N
    idx = int(rand() * nstars) + 1

    do j = idx, nstars - 1
        
        x(j) = x(j + 1)
        y(j) = y(j + 1)
        z(j) = z(j + 1)

        vx(j) = vx(j + 1)
        vy(j) = vy(J + 1)
        vz(j) = vz(j + 1)
    
    end do

    nstars = nstars - 1
end do

contains
!
!random number generator
!
real function rand()
real :: rand_val
call random_number(rand_val)
rand = rand_val
end function rand


end subroutine pruning_excess