program poisson_mpi
	use poisson_module
	use mpi
	implicit none

	integer :: i,j,k,ret
	integer :: id1,id2,id3,id4,id5
	real, allocatable, dimension(:) :: u,u0,f
	integer, allocatable, dimension(:) :: buf
	real :: h, l2
	integer :: prank, psize, ierror, tag, status(MPI_STATUS_SIZE)
	double precision :: t_start, t_end
	integer :: n_loc, offset;
	integer :: i_start, i_end

! Initialization of the MPI library
	call MPI_INIT(ierror)	
	call MPI_COMM_RANK(MPI_COMM_WORLD, prank, ierror)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, ierror)


	h = 1. / (N+1)
! Divide n in psize and distribute the excess to the (n % psize) proc
	if (prank < mod((N+1),psize)) then
		n_loc = ((N+1) / psize) + 1
	else
		n_loc = ((N+1) / psize)
	endif
! Computing the offset of where in the global array the local array is
! located. (this is needed to initialize f) It could be computed locally
! without communication but this way you see a Allgather
	allocate(buf(psize))
	call MPI_ALLGATHER(n_loc, 1, MPI_INT, buf, 1, MPI_INT, MPI_COMM_WORLD,ierror)
	offset = 0
	do i = 1, prank
		offset = offset + buf(i)
	enddo
! add 2 for north and south ghost
	n_loc = n_loc + 2

	allocate(u((N+1)*n_loc),u0((N+1)*n_loc),f((N+1)*n_loc))

	u = 0.
	u0 = 0.
	do i = 2, N
		do j = 2, n_loc
			id1 = id(n_loc,i,j)
			f(id1) = -2.*100. * M_PI * M_PI * sin(10.*M_PI*(i-1)*h) * sin(10.*M_PI*(j-1)*h)
		enddo
	enddo	


	k = 0
	l2 = 1.
	t_start = MPI_WTIME()

!	do while (l2 > eps)
	do while (k < 10)
		l2 = 0.
		i_start = 0
		i_end = 0

!    if(prank > 0) { // send recv with top proc
!      MPI_Sendrecv(uo[1], n, MPI_FLOAT, prank - 1, NORTH,
!                   uo[0], n, MPI_FLOAT, prank - 1, SOUTH, MPI_COMM_WORLD, &status);
!    } else {
!      ++i_start;
!    }
!
!    if(prank < psize - 1) { // send recv with bottom proc
!      MPI_Sendrecv(uo[n_loc - 2], n, MPI_FLOAT, prank + 1, SOUTH,
!                   uo[n_loc - 1], n, MPI_FLOAT, prank + 1, NORTH, MPI_COMM_WORLD, &status);
!    } else {
!      --i_end;
!    }

! send recv with top proc
		if (prank > 0)  then 
			call MPI_SENDRECV(u0(1:N+1), (N+1), MPI_FLOAT, prank - 1, NORTH, u0(1:N+1), (N+1), MPI_FLOAT, prank - 1, SOUTH, MPI_COMM_WORLD, status, ierror)
		else
			i_start = i_start + 1
		endif

! send recv with bottom proc
		if(prank < psize - 1) then
			call MPI_SENDRECV(u0(n_loc - 2), (N+1), MPI_FLOAT, prank + 1, SOUTH, u0(n_loc - 1), (N+1), MPI_FLOAT, prank + 1, NORTH, MPI_COMM_WORLD, status, ierror)
		else
			i_end = i_end - 1
		endif

		do i = i_start, i_end
			do j = 2, N
				id1 = id(N,i,j)
				id2 = id(N,i-1,j)
				id3 = id(N,i+1,j)
				id4 = id(N,i,j-1)
				id5 = id(N,i,j+1)
				u(id1) = 0.25 * ( u0(id2) + u0(id3) + u0(id4) + u0(id5) - f(id1)*h*h);
				! L2 norm
				l2 = l2 + (u0(id1) - u(id1))*(u0(id1) - u(id1));
			enddo
		enddo

		if (prank.eq.0) call print_mat(u,n_loc,(N+1))

		u0 = u

! reduce the 12 norm to every proc
		call MPI_ALLREDUCE(MPI_IN_PLACE, l2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierror);


		! output	
		if (prank .eq. 0) then
			write(*,*) 'iteration',k,sqrt(l2),prank,i_start,i_end,'size of matrix',(N+1),n_loc
		endif
!			ret=write_to_file_binary(N,c_loc(u(1)),k,-1.,1.)
		k = k + 1
	enddo	
	t_end = MPI_WTIME()
	if (prank .eq. 0) then
		write(*,*) 'T',(t_end - t_start),'[s]'
	endif

	deallocate(u,u0,f)

	call MPI_FINALIZE(ierror)

end program poisson_mpi

