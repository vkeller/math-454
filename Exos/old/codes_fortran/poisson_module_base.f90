module poisson_module

	use ISO_C_BINDING
	implicit none

	integer, parameter :: NORTH = 0, SOUTH = 1
	integer, parameter :: N=15
	real, parameter :: eps=.005
	real, parameter :: M_PI  = 2.D0*DASIN(1.D0)

public :: id

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

contains


!* --------------------------------------
!* Write images with a 1D array
!* --------------------------------------


	subroutine write_to_file_binary1D(siz, tmptab, iter, minval, maxval)
		implicit none
		
		integer, intent(in) :: siz
		real, dimension(:), intent(in) :: tmptab
		integer, intent(in) :: iter	
		real :: minval, maxval
! local
		integer :: w,h
		character*12  :: fnameout
		integer :: i, j, k,x,y
		integer :: itmp, icnt
		character*14 :: frmtstr
		character*54 :: headmsw
		character*4 :: byt4
		character*2 :: byt2
		character, dimension(:), allocatable :: img
		real :: v,r,g,b
		integer ::  padding

		w = siz
		h = siz

		itmp = mod(w, 4)
!		if (itmp .NE. 0) then
!			write(*,*) 'width must be multiple of 4'
!			stop
!		endif

		allocate(img(3*w*h))
		padding = mod((4 - mod((w*3),4)),4);

		do i = 1, w
			do j = 1, h
!				x = (i-1)
!				y = (w - 1) - (j-1)
				x = i
				y = w - j

				v = ((tmptab(w*(j-1)+i) - minval)/(maxval - minval))
				r = v * 255 ! Red channel
				g = v * 255 ! Green channel
				b = v * 255 ! Blue channel

				r = min(r, 255.)
				g = min(g, 255.)
				b = min(b, 255.)

				img((x + y*w)*3 + y*padding + 2) = char(int(r));
				img((x + y*w)*3 + y*padding + 1) = char(int(g));
				img((x + y*w)*3 + y*padding + 0) = char(int(b));
			enddo
		enddo

		write(fnameout,'(''out_'',i3.3,''.bmp'')') iter ! name of BMP file
		open(unit=2,file=fnameout,status='unknown')
!* header 1 (file header ; 1--14 byte)
		headmsw( 1: 2) = 'BM'             ! declaring this is BMP file
		itmp = 54 + (w * h * 3) ! total file size = header + data
		call num2bit4(itmp,byt4)
		headmsw( 3: 6) = byt4(1:4)
		itmp = 0                        ! may be 0
		call num2bit2(itmp,byt2)
		headmsw( 7: 8) = byt2(1:2)
		itmp = 0                        ! may be 0
		call num2bit2(itmp,byt2)
		headmsw( 9:10) = byt2(1:2)
		itmp = 54                       ! must be 54 : total length of header
		call num2bit4(itmp,byt4)
		headmsw(11:14) = byt4(1:4)
!* header 2 (bit-map header ; 13--54 byte)
		itmp = 40                       ! must be 40 : length of bit-map header
		call num2bit4(itmp,byt4)
		headmsw(15:18) = byt4(1:4)
		itmp = w          	         ! width
		call num2bit4(itmp,byt4)
		headmsw(19:22) = byt4(1:4)
		itmp = h	                   ! height
		call num2bit4(itmp,byt4)
		headmsw(23:26) = byt4(1:4)
		itmp = 1                        ! must be 1
		call num2bit2(itmp,byt2)
		headmsw(27:28) = byt2(1:2)
		itmp = 24                       ! must be 24 : color depth in bit.
		call num2bit2(itmp,byt2)
		headmsw(29:30) = byt2(1:2)
		itmp = 0                        ! may be 0 : compression method index
		call num2bit4(itmp,byt4)
		headmsw(31:34) = byt4(1:4)
		itmp = 0                        ! may be 0 : file size if compressed
		call num2bit4(itmp,byt4)
		headmsw(35:38) = byt4(1:4)
		itmp = 0                        ! arbit. : pixel per meter, horizontal
		call num2bit4(itmp,byt4)
		headmsw(39:42) = byt4(1:4)
		itmp = 0                        ! arbit. : pixel per meter, vertical
		call num2bit4(itmp,byt4)
		headmsw(43:46) = byt4(1:4)
		itmp = 0                        ! may be 0 here : num. of color used
		call num2bit4(itmp,byt4)
		headmsw(47:50) = byt4(1:4)
		itmp = 0                        ! may be 0 here : num. of important color
		call num2bit4(itmp,byt4)
		headmsw(51:54) = byt4(1:4)

!* writing header part
		write(2,'(a54,$)') headmsw(1:54)
!* image data
		itmp = (w * h * 3) 
		write(frmtstr,'(''('',i8.8,''A,$)'')') itmp
		write(2,fmt=frmtstr) img
		close(2)
		deallocate(img)

		return
	end subroutine write_to_file_binary1D

!* --------------------------------------
!* Write images with a 2D array
!* --------------------------------------

	subroutine write_to_file_binary2D(siz, tmptab, iter, minval, maxval)
		implicit none
		
		integer, intent(in) :: siz
		real, dimension(:,:), intent(in) :: tmptab
		integer, intent(in) :: iter	
		real :: minval, maxval
! local
		integer :: w,h
		character*12  :: fnameout
		integer :: i, j, k,x,y
		integer :: itmp, icnt
		character*14 :: frmtstr
		character*54 :: headmsw
		character*4 :: byt4
		character*2 :: byt2
		character, dimension(:), allocatable :: img
		real :: v,r,g,b
		integer ::  padding

		w = siz
		h = siz

		itmp = mod(w, 4)
!		if (itmp .NE. 0) then
!			write(*,*) 'width must be multiple of 4'
!			stop
!		endif

		allocate(img(3*w*h))
		padding = mod((4 - mod((w*3),4)),4);

		do i = 1, w
			do j = 1, h
!				x = (i-1)
!				y = (w - 1) - (j-1)
				x = i
				y = w - j

				v = ((tmptab(j,i) - minval)/(maxval - minval))
				r = v * 255 ! Red channel
				g = v * 255 ! Green channel
				b = v * 255 ! Blue channel

				r = min(r, 255.)
				g = min(g, 255.)
				b = min(b, 255.)

				img((x + y*w)*3 + y*padding + 2) = char(int(r));
				img((x + y*w)*3 + y*padding + 1) = char(int(g));
				img((x + y*w)*3 + y*padding + 0) = char(int(b));
			enddo
		enddo

		write(fnameout,'(''out_'',i3.3,''.bmp'')') iter ! name of BMP file
		open(unit=2,file=fnameout,status='unknown')
!* header 1 (file header ; 1--14 byte)
		headmsw( 1: 2) = 'BM'             ! declaring this is BMP file
		itmp = 54 + (w * h * 3) ! total file size = header + data
		call num2bit4(itmp,byt4)
		headmsw( 3: 6) = byt4(1:4)
		itmp = 0                        ! may be 0
		call num2bit2(itmp,byt2)
		headmsw( 7: 8) = byt2(1:2)
		itmp = 0                        ! may be 0
		call num2bit2(itmp,byt2)
		headmsw( 9:10) = byt2(1:2)
		itmp = 54                       ! must be 54 : total length of header
		call num2bit4(itmp,byt4)
		headmsw(11:14) = byt4(1:4)
!* header 2 (bit-map header ; 13--54 byte)
		itmp = 40                       ! must be 40 : length of bit-map header
		call num2bit4(itmp,byt4)
		headmsw(15:18) = byt4(1:4)
		itmp = w          	         ! width
		call num2bit4(itmp,byt4)
		headmsw(19:22) = byt4(1:4)
		itmp = h	                   ! height
		call num2bit4(itmp,byt4)
		headmsw(23:26) = byt4(1:4)
		itmp = 1                        ! must be 1
		call num2bit2(itmp,byt2)
		headmsw(27:28) = byt2(1:2)
		itmp = 24                       ! must be 24 : color depth in bit.
		call num2bit2(itmp,byt2)
		headmsw(29:30) = byt2(1:2)
		itmp = 0                        ! may be 0 : compression method index
		call num2bit4(itmp,byt4)
		headmsw(31:34) = byt4(1:4)
		itmp = 0                        ! may be 0 : file size if compressed
		call num2bit4(itmp,byt4)
		headmsw(35:38) = byt4(1:4)
		itmp = 0                        ! arbit. : pixel per meter, horizontal
		call num2bit4(itmp,byt4)
		headmsw(39:42) = byt4(1:4)
		itmp = 0                        ! arbit. : pixel per meter, vertical
		call num2bit4(itmp,byt4)
		headmsw(43:46) = byt4(1:4)
		itmp = 0                        ! may be 0 here : num. of color used
		call num2bit4(itmp,byt4)
		headmsw(47:50) = byt4(1:4)
		itmp = 0                        ! may be 0 here : num. of important color
		call num2bit4(itmp,byt4)
		headmsw(51:54) = byt4(1:4)

!* writing header part
		write(2,'(a54,$)') headmsw(1:54)
!* image data
		itmp = (w * h * 3) 
		write(frmtstr,'(''('',i8.8,''A,$)'')') itmp
		write(2,fmt=frmtstr) img
		close(2)
		deallocate(img)

		return
	end subroutine write_to_file_binary2D

!* --------------------------------------
!* convert integer values to 4 8-bit characters
!* --------------------------------------

	subroutine num2bit4(inum,byt4)
		implicit none
		integer inum
		character*4 byt4
		integer itmp1, itmp2
		itmp1 = inum
		itmp2 = itmp1 / 256**3
		byt4(4:4) = char(itmp2)
		itmp1 =-itmp2 * 256**3 +itmp1
		itmp2 = itmp1 / 256**2
		byt4(3:3) = char(itmp2)
		itmp1 =-itmp2 * 256**2 +itmp1
		itmp2 = itmp1 / 256
		byt4(2:2) = char(itmp2)
		itmp1 =-itmp2 * 256    +itmp1
		byt4(1:1) = char(itmp1)
		return
	end subroutine num2bit4

!* --------------------------------------
!* convert integer values to 2 8-bit characters
!* --------------------------------------

	subroutine num2bit2(inum,byt2)
		implicit none
		integer inum
		character*2 byt2
		integer itmp1, itmp2
		itmp1 = inum
		itmp2 = itmp1 / 256
		byt2(2:2) = char(itmp2)
		itmp1 =-itmp2 * 256 + itmp1
		byt2(1:1) = char(itmp1)
		return
	end subroutine num2bit2


end module poisson_module
