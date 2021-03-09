program poisson 
	use poisson_module
	
	implicit none

	integer, parameter :: N=256
	real, parameter :: eps=.005
	real, parameter :: M_PI  = 2.D0*DASIN(1.D0)


	integer :: i,j,k,ret
	real, allocatable, dimension(:,:) :: u,u0,f
	real :: h, l2

	h = 1. / (N+1)

	allocate(u(N+1,N+1),u0(N+1,N+1),f(N+1,N+1))
	
	u0 = 0.
	do i = 2, N
		do j = 2, N
			f(i,j) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*(i-1)*h) * sin(10.*M_PI*(j-1)*h)
		enddo
	enddo

	k = 0
	l2 = 1.
	do while (l2 > eps)
		l2 = 0.

		do i = 2, N
			do j = 2, N
				u(i,j) = 0.25 * ( u0(i-1,j) + u0(i+1,j) + u0(i,j-1) + u0(i,j+1) - f(i,j)*h*h);
				! L2 norm
				l2 = l2 + (u0(i,j) - u(i,j))*(u0(i,j) - u(i,j));
			enddo
		enddo

		u0 = u
		! output	
		write(*,*) 'iteration',k,sqrt(l2)
		ret=write_to_file_binary(N,c_loc(u(1,1)),k,-1.,1.)
		k = k + 1
	enddo	
	deallocate(u,u0,f)

end program poisson
