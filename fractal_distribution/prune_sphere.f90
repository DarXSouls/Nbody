subroutine pruning_sphere(x, y, z, vx, vy, vz, nstars, N, rad_cluster)

integer, intent(inout) :: nstars
integer, intent(inout) :: N

real, dimension(:), intent(inout) :: x, y, z
real, dimension(:), intent(inout) :: vx, vy, vz

integer :: i, j

double precision :: rad_cluster
real :: dist

integer, dimension(:), allocatable :: to_remove
integer :: remove_count

allocate(to_remove(nstars))

remove_count = 0

do i = 1, nstars

  dist = sqrt((x(i) - 0.5)**2 + (y(i) - 0.5)**2 + (z(i) - 0.5)**2)

    if (dist > rad_cluster) then
      remove_count = remove_count + 1
      to_remove(remove_count) = 1
    end if
end do

do i = remove_count, 1, -1

    do j = 1, nstars - 1
         
        x(j) = x(j + 1)
        y(j) = y(J + 1)
        z(j) = z(j + 1)

        vx(j) = vx(j + 1)
        vy(j) = vy(j + 1)
        vz(j) = vz(j + 1)

    end do

    nstars = nstars - 1

end do

deallocate(to_remove    )

end subroutine pruning_sphere