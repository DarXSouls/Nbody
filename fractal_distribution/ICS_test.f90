program gen_ICS

use global_variables
use global_parameters

implicit none

character*100 :: filename
integer filename_length
integer :: i


call particle_generation


print *, 'filename ?'
read(*,*) filename

filename_length = len_trim(filename)


open(71, file= "ics.dat", status='old')

allocate(x(N), y(N), z(N))
allocate(vx(N), vy(N), vz(N))
allocate(m(N))

do i = 1, N
    read(71, *) x(i), y(i), z(i), vx(i), vy(i), vz(i), m(i) 
end do

close(71)




open(70, file=filename(1:filename_length), form='unformatted')

global_time = 0

! write data

write(70) N
write(70) global_time
write(70) units_mass_g, units_length_cm, units_time_s

write(70) x(1:N)
write(70) y(1:N)
write(70) z(1:N)
write(70) vx(1:N)
write(70) vy(1:N)
write(70) vz(1:N)
write(70) m(1:N)

close(70)


deallocate(x, y, z, vx, vy, vz, m)


end program gen_ICS