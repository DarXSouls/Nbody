module timing_info

!
! coding timing variables
!
integer*8 timing_cr, timing_cm
double precision :: timing_rate
integer*8 :: timing_start, timing_end ! start and end of code
integer*8 :: t0, t1 ! generic start/end times
integer*8 :: th0, th1 ! generic start/end times
integer*8 :: tw0, tw1 ! generic start/end times
! various parts we want to time
double precision :: timing_accel, timing_dts, timing_hermite, timing_herm2, time_so_far 
double precision :: timing_snaps, timing_diag
integer*8 :: accel_function_counter ! counts the total number of pair-wise force calculations
integer*8 :: calls_to_get_accel     ! counts the total number of calls to get_accel_*() 
end module timing_info
