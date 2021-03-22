program LU_decomposition_TDM
	implicit none
	double precision, allocatable, dimension(:) :: xnode,f
	double precision, allocatable, dimension(:) :: a,b,c,L,U,y,z,x,dfdx_analytical
	integer :: n,i
	double precision :: h
	
	n=26 	!no. of grid points
	allocate(xnode(1:n),f(1:n))
	allocate(a(1:n),b(1:n),c(1:n),L(1:n),U(1:n))
	allocate(y(1:n),z(1:n),x(1:n),dfdx_analytical(1:n))

	h=3.0d0/(n-1)
	
	!xgrid values:
	xnode(1) = 0.0d0
	do i=2,n,1
		xnode(i) = xnode(i-1) + h
	end do
	!Function values:
	do i=1,n,1
		f(i) = sin(5.0d0*xnode(i))
	end do
	
	!RHS vector:
	do i=1,n,1
		if(i==1) then
			y(i) = (1.0d0/h)*((-5.0d0/2.0d0)*f(i) + 2.0d0*f(i+1) + (1.0d0/2.0d0)*f(i+2))
		elseif(i==n) then
			y(i) = (1.0d0/h)*((5.0d0/2.0d0)*f(i) - 2.0d0*f(i-1) - (1.0d0/2.0d0)*f(i-2))
		else
			y(i) = (3.0d0/h)*(f(i+1)-f(i-1))
		end if
	end do
	
	!digonal elements
	b(:) = 4.0d0
	b(1) = 1.0d0
	b(n) = 1.0d0
	
	!sub-diagonal elements
	a(:) = 1.0d0
	a(1) = 0.0d0
	a(n) = 2.0d0
	
	!super-diagonal elements
	c(:) = 1.0d0
	c(1) = 2.0d0
	c(n) = 0.0d0
	
	!Upper triangular and Lower triangular elements:
	U(1) = b(1)
	do i = 1,n,1
		L(i+1) = a(i+1)/U(i)
		U(i+1) = b(i+1) - L(i+1)*c(i)
	end do
	
	!Forward substitution
	z(1) = y(1)
	do i = 2,n,1
		z(i) = y(i) - L(i)*z(i-1)
	end do
	
	!Backward substitution
	x(n) = z(n)/U(n)
	do i = n-1,1,-1
		x(i) = (z(i) - c(i)*x(i+1))/U(i)
	end do
	
	!Analytical solution:
	do i=1,n,1
		dfdx_analytical(i) = 5.0d0*cos(5.0d0*xnode(i))
	end do
	
	!Printing and writing the values:
	do i=1,n,1
		write(*,*),xnode(i),x(i),dfdx_analytical(i)
	end do
	
	open(unit=10, file='first_derivative.txt')
		do i=1,n,1
			write(10,*),xnode(i),x(i),dfdx_analytical(i)
		end do
	close(10)


end program LU_decomposition_TDM
