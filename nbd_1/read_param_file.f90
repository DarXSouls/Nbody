subroutine read_param_file

!
!
! Reads the parameter file. Note that 'parameters' here are actually fortran
! variables, since they are not set at compile time (unlike true fortran parameters).
!
! parmeters in file must have the following format:
!   parameter_name
!   parameter_value
!

use global_variables
use global_parameters

character*30 parameter_name
integer :: param_name_len

!
! set defaults
!

system_dt = -1 ! i.e. force system to evolve on single dt to be off by default


!
! open the file. For now, the file always has the same name, but can change this
! in the future by using command-line arguments.
!
open(10, file='nbd.param')

!
! read lines of parameter file
!
do
   !
   ! read parameter name; break if EoF
   !
   read(10, *, IOSTAT=io) parameter_name
   if (io.ne.0) exit
   !
   ! different cases depending on parameter name
   !
   param_name_len = len_trim(parameter_name)
   select case( parameter_name(1:param_name_len) )
      case('initial_conditions')
         read(10, *, IOSTAT=io) initial_conditions
      case('end_time')
         read(10, *, IOSTAT=io) end_time
      case('dt_max')
         read(10, *, IOSTAT=io) dt_max
      case('int_mode')
         read(10, *, IOSTAT=io) int_mode
      case('snapshot_base')
         read(10, *, IOSTAT=io) snapshot_base
      case('time_between_snapshots')
         read(10, *, IOSTAT=io) time_between_snapshots
      case('softening')
         read(10, *, IOSTAT=io) soft
      case('restart_flag')
         read(10, *, IOSTAT=io) restart_flag
      case('diagnostic_freq')
         read(10, *, IOSTAT=io) diagnostic_freq
      case('system_dt')
         read(10, *) system_dt
      case default
         print *, "Unrecognised parameter ", parameter_name(1:param_name_len)
         print *, "Please check and try again!"
         stop
   end select
   if (io.ne.0) then
      print *, "missing/incorrect value for parameter ", parameter_name(1:param_name_len)
      stop
   end if
end do

!
! check consistency of options
!
if (system_dt .gt. 0  .and. int_mode .ne. 2) then
   print *, "system_dt option is only valid with block dts (int_mode = 2)!"
   print *, "stopping"
   stop
end if

!
! close the parameter file
!
close(10)

end subroutine read_param_file
