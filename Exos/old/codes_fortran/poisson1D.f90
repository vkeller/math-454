program poisson1D
	use poisson_module
	
	implicit none

	integer :: i,j,k,ret,n_min
	real, allocatable, dimension(:) :: u,u0,f
	real :: h, l2
	real :: t_start, t_end

	n_min = (N+1)
	h = 1. / n_min

	allocate(u(n_min*n_min),u0(n_min*n_min),f(n_min*n_min))
	
	u0 = 0.
	do i = 2, N
		do j = 2, N
			f(N*i+j) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*(i-1)*h) * sin(10.*M_PI*(j-1)*h)
		enddo
	enddo

	k = 0
	l2 = 1.
	call CPU_TIME(t_start)
	do while (l2 > eps)
		l2 = 0.
		! width*row+col
		do i = 2, N
			do j = 2, N
				u(N*i+j) = 0.25 * ( u0(N*(i-1)+j) + u0(N*(i+1)+j) + u0(N*i+(j-1)) + u0(N*i+(j+1)) - f(N*i+j)*h*h);
				! L2 norm
				l2 = l2 + (u0(N*i+j) - u(N*i+j))*(u0(N*i+j) - u(N*i+j));
			enddo
		enddo

		u0 = u
		! output	
!		write(*,*) 'iteration',k,sqrt(l2)
!		ret=write_to_file_binary(N,c_loc(u(1)),k,-1.,1.)
		k = k + 1
	enddo	
	call CPU_TIME(t_end)
	write(*,*) 'T',(t_end - t_start),'[s]'
	write(*,*) 'Nb steps',k,'l2',l2
	deallocate(u,u0,f)

end program poisson1D



