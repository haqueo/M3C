module sync
	implicit none
	complex(kind=16), parameter :: ii=cmplx(0.0,1.0) !ii = sqrt(-1)
	integer :: ntotal !total number of nodes
	real(kind=8) :: c !coupling coefficient
        real(kind=8), allocatable, dimension(:) :: w !array of frequencies
	save
    contains


!---------------------------
subroutine rk4(t0,y0,dt,nt,y,order)
    !4th order RK method
    !input:
    !t0: initial time
    !y0: initial condition (array)
    !dt: time step
    !nt: number of time steps
    !output:
    !y: solution at each time step (array)
    !order: order parameter
    implicit none
    real(kind=8), dimension(:), intent(in) :: y0
    real(kind=8), intent(in) :: t0,dt
    integer, intent (in) :: nt
    real(kind=8), dimension(size(y0),nt+1), intent(out) :: y
    real(kind=8), dimension(nt), intent(out) :: order
    real(kind=8), dimension(size(y0)) :: f1, f2, f3, f4
    real(kind=8) :: t,halfdt,fac,sumCosines,sumSines
    integer:: k,N,i1
    
        halfdt = 0.5d0*dt
        fac = 1.d0/6.d0

        y(:,1) = y0 !initial condition
        t = t0 !initial time

        N = size(y0)
        do k = 1, nt !advance nt time steps

           f1 = dt*RHS(t, y(:,k))

           f2 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f1)

           f3 = dt*RHS(t + halfdt, y(:,k) + 0.5d0*f2)

           f4 = dt*RHS(t + dt, y(:,k) + f3)

           y(:,k+1) = y(:,k) + (f1 + 2*f2  + 2*f3 + f4)*fac

           t = t + dt
        
        !add code to compute order
        
        sumSines = dble(0)
        sumCosines = dble(0)
        
        DO i1 = 1,N
            sumCosines = sumCosines + cos(y(i1,k))
            sumSines = sumSines + sin(y(i1,k))
        END DO
        
        order(k) = dble(SQRT((sumCosines**2)+(sumSines**2)))/dble(N)
           
        end do
        
        
end subroutine rk4

!---------------------------
function RHS(t,f)
    !RHS sync
    !f is the array of phases at time, t
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    real(kind=8), dimension(size(f)) :: RHS
    integer :: i1,N,j
    real(kind=8) :: sumsines
    
    N = size(f)
   
    !w, is a module variable, get the normals from python
    !c is also a module variable
    DO i1 = 1, N
        sumsines = dble(0)
        DO j = 1,N
            sumsines = sumsines + sin(f(i1)-f(j))
        END DO
    
    
        
        RHS(i1) = w(i1) - (c/N)*sumsines
    END DO
 
 
end function RHS
!---------------------------
end module sync


