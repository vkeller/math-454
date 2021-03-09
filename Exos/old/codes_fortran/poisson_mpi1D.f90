program poisson_mpi1D
	use poisson_module
	use mpi
	implicit none

	integer :: i,j,k,ret
	real, allocatable, dimension(:) :: u,u0,f,u_global
	integer, allocatable, dimension(:) :: buf
	real :: h, l2
	integer :: prank, psize, ierror, tag, status(MPI_STATUS_SIZE)
	double precision :: t_start, t_end
	integer :: n_loc, offset, n_min;
	integer :: j_start, j_end, i_start, i_end
	real :: x,y

! Initialization of the MPI library
	call MPI_INIT(ierror)	
	call MPI_COMM_RANK(MPI_COMM_WORLD, prank, ierror)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, ierror)

	n_min = (N+1)
	h = 1. / n_min
! Divide n in psize and distribute the excess to the (n % psize) proc
	if (prank < mod(n_min,psize)) then
		n_loc = (n_min / psize) + 1
	else
		n_loc = (n_min / psize)
	endif
! Computing the offset of where in the global array the local array is
! located. (this is needed to initialize f) It could be computed locally
! without communication but this way you see a Allgather
	allocate(buf(psize))
	call MPI_ALLGATHER(n_loc, 1, MPI_INTEGER, buf, 1, MPI_INTEGER, MPI_COMM_WORLD,ierror)
	offset = 0
	do i = 1, prank
		offset = offset + buf(i)
	enddo
! add 2 for east and west ghost
	n_loc = n_loc + 2

	allocate(u(n_min*n_loc),u0(n_min*n_loc),f(n_min*n_loc))
	if (prank.eq.0) allocate(u_global(n_min*n_min))

	u = 0.
	u0 = 0.
	f = 0.


	do i = 2, n_min-1
		do j = 2, n_loc-1
!			x = real(i) * h
!			y = real(j - j_start + offset) * h
!			f(n_loc*i+j) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*x) * sin(10.*M_PI*y)
			f(n_loc*i+j) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*(i-1)*h) * sin(10.*M_PI*(j-1)*h)

		enddo
	enddo	
	
!if (prank.eq.1) write(*,*) prank, f

	k = 0
	l2 = 1.
	t_start = MPI_WTIME()


!	do while (l2 > eps)
	do while (k < 400)
		l2 = 0.
		j_start = 2
		j_end   = n_loc-1
		i_start = 2
		i_end   = n_min-1

! send recv with left proc
		if (prank > 0)  then 
			call MPI_SENDRECV(u0(n_min+2), n_min, MPI_REAL, prank - 1, NORTH, u0(1), n_min, MPI_REAL, prank - 1, SOUTH, MPI_COMM_WORLD, status, ierror)
!			call MPI_SENDRECV(u0(n_min:2*n_min), n_min, MPI_REAL, prank - 1, NORTH, u0(1:n_min), n_min, MPI_REAL, prank - 1, SOUTH, MPI_COMM_WORLD, status, ierror)
!				call MPI_SENDRECV(u0(:,2), n_min, MPI_REAL, prank - 1, NORTH, u0(:,1), n_min, MPI_REAL, prank - 1, SOUTH, MPI_COMM_WORLD, status, ierror)
		endif

! send recv with right proc
		if(prank < psize - 1) then
			call MPI_SENDRECV(u0(n_loc-1), n_min, MPI_REAL, prank + 1, SOUTH, u0(n_loc), n_min, MPI_REAL, prank + 1, NORTH, MPI_COMM_WORLD, status, ierror)
!			call MPI_SENDRECV(u0((n_loc-1)*n_min:(n_loc-1)*2*n_min), n_min, MPI_REAL, prank + 1, SOUTH, u0(n_loc:n_loc+n_min), n_min, MPI_REAL, prank + 1, NORTH, MPI_COMM_WORLD, status, ierror)
!			call MPI_SENDRECV(u0(:,n_loc-1), n_min, MPI_REAL, prank + 1, SOUTH, u0(:,n_loc), n_min, MPI_REAL, prank + 1, NORTH, MPI_COMM_WORLD, status, ierror)
		endif

		do i = i_start, i_end
			do j = j_start,j_end
				u(n_loc*i+j) = 0.25 * ( u0(n_loc*(i-1)+j) + u0(n_loc*(i+1)+j) + u0(n_loc*i+(j-1)) + u0(n_loc*i+(j+1)) - f(n_loc*i+j)*h*h);
				! L2 norm
				l2 = l2 + (u0(n_loc*i+j) - u(n_loc*i+j))*(u0(n_loc*i+j) - u(n_loc*i+j));
			enddo
		enddo


! reduce the 12 norm to every proc
		call MPI_ALLREDUCE(MPI_IN_PLACE, l2, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierror);

		call MPI_Gather(u(1),(n_min*(n_loc-2)),MPI_REAL, u_global,(n_min*(n_loc-2)) ,  MPI_REAL, 0, MPI_COMM_WORLD,ierror)
!		if (prank.eq.0) call write_to_file_binary1D(n_min, u_global, k, -1., 1.)

		! output
		if (prank .eq. 0) then
			write(*,*) 'iteration',k,l2
		endif
		k = k + 1
		u0 = u
	enddo	
	t_end = MPI_WTIME()
	if (prank .eq. 0) then
		write(*,*) 'T',(t_end - t_start),'[s]'
		write(*,*) 'Nb steps',k,'l2',l2
		write(*,*) 'Size',n_min
	endif

	deallocate(u,u0,f)
	if (prank.eq.0) deallocate(u_global)

	call MPI_FINALIZE(ierror)

end program poisson_mpi1D
