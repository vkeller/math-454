module poisson_module_c

	use ISO_C_BINDING
	implicit none

! =====================================================================================
! ISO_C_BINDING interface
!
! C function prototype :
! int write_to_file_binary(int n, float **tmptab, int iter, float minval, float maxval);
! =====================================================================================
interface
	integer (C_INT) function write_to_file_binary &
		(n, tmptab, iter, minval, maxval) &
		bind(C, NAME='write_to_file_binary1D')
		use ISO_C_BINDING
		implicit none
		integer(c_int),value :: n
		type (c_ptr), value :: tmptab
		integer(c_int) ,value :: iter
		real(c_float), value :: minval
		real(c_float), value :: maxval
	end function write_to_file_binary
end interface

end module poisson_module_c
