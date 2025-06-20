subroutine get_ics

use global_variables
use global_parameters

integer :: io_stat

!
! Open the file and read in the data. The snapshots use the same format as the ics,
! ics, excepts that they hold additional information in a footer
!
open(20, file = initial_conditions, form='unformatted', status='old', iostat=io_stat)

if(io_stat /= 0) then
    print *, 'Error opening the file:', initial_conditions
    stop
end if


!
! First read in the number of bodies and other header info
!
read(20, iostat=io_stat) N
print *, 'read N!', N

if(io_stat /= 0) then
   print *, 'Error reading N'
   stop
end if

read(20) global_time
print *, 'read global time ', global_time

read(20) units_mass_g, units_length_cm, units_time_s
print *, 'Units of mass_g: ', units_mass_g
print *, 'Units of length: ', units_length_cm
print *, 'Units of time: ', units_time_s

!
! now allocate all the arrays for the program
!
call allocate_memory

!
! read in the arrays
!
read(20, iostat=io_stat) x(1:N)

if (io_stat /= 0) then
   print *, 'Error readin x from file'
   stop
end if 

read(20) y(1:N)
read(20) z(1:N)
read(20) vx(1:N)
read(20) vy(1:N)
read(20) vz(1:N)
read(20) m(1:N)
if (system_dt .gt. 0) then
read(20) is_system(1:N)
end if
!
! close the file
!
close(20)

print *, 'pos 1', x(1), y(1), z(1)
print *, 'pos 2', x(2), y(2), z(2)
print *, 'vel 1', vx(1), vy(1), vz(1)
print *, 'vel 2', vx(2), vy(2), vz(2)

end subroutine get_ics
