program poisson2D
	use poisson_module
	
	implicit none

	integer :: i,j,k,ret
	real, allocatable, dimension(:,:) :: u,u0,f
	real :: h, l2
	real :: t_start, t_end

	h = 1. / (N+1)

	allocate(u((N+1),(N+1)),u0((N+1),(N+1)),f((N+1),(N+1)))
	
	u0 = 0.
	do i = 2, N
		do j = 2, N
			f(i,j) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*(i-1)*h) * sin(10.*M_PI*(j-1)*h)
		enddo
	enddo

	k = 0
	l2 = 1.
	call CPU_TIME(t_start)
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
!		write(*,*) 'iteration',k,sqrt(l2)
!		ret=write_to_file_binary(N,c_loc(u(1)),k,-1.,1.)
		k = k + 1
	enddo	
	call CPU_TIME(t_end)
	write(*,*) 'T',(t_end - t_start),'[s]'
	write(*,*) 'Nb steps',k,'l2',l2

	deallocate(u,u0,f)

end program poisson2D
