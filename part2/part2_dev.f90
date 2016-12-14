module flunet
        use omp_lib
	implicit none
        !add variables as needed
	save
	contains


subroutine rhs(n,y,t,a,b0,b1,g,k,w,dy,qnet,A)
    implicit none
    !Return RHS of network flu model
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !output: dy, RHS
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n) :: S,E,C
    integer :: i
    real(kind=8), dimension(n),intent(in) :: qnet
    real(kind=8), dimension(n,n),intent(in) :: A
    
    S = y(1:n)
    E = y(n+1:2*n)
    C = y((2*n)+1:3N)
    
    b = b0 + b1*(1.d0+cos(2.0*.pi*t))
    
    !need use P values without creating an NxN matrix
    
    
    !Choices: (i) add P as an input variable and use it
    !         (ii) add qnet, A as input variables then use them
    ! Considering we need qnet, A to make P anyway - (ii) is probably cheaper
    !but this still uses an NxN matrix..
    !yeah, this is cheaper since we can use vectorisation.
    
    
    
    !from python:
    !P[i,j] = float(qnet[i]*A[i,j])/float(np.dot(qnet,A[:,j]))
    !fortran equivalent:
    !P(i,j) = dble(qnet(i)*A(i,j))/float(dot_product(qnet,A(:,j)))
    
    
    !Now, what we want is terms of the form: P(:,i)*dotproduct*S
    !P(:,i)*DOTPRODUCT*S = [qnet(1)*A(1,i)/dot_product(qnet,A(:,i)),...,qnet(N)*A(N,i)/dot_product(qnet,A(:,N))]*dotproduct*S
    !                    = (1/dot_product(qnet,A(:,i))[qnet(1)*A(1,i),qnet(2)*A(1,i),..,qnet(N)*A(N,i)]*dotproduct*S
    !                    = (1/dot_product(qnet,A(:,i))*(qnet*elementwisemultiplication*A(:,i))*dotproduct*S
    !similar for E,C.
    
    
    
    
    DO i = 1:n
        !dSi/dt 
        dy(i) = k*(1-S(i))-(b*C(i)*S(i))-(w*S(i))+w*(1.d0/dble(dot_product(qnet,A(:,i))))*(dot_product(qnet*A(:,i),S)) !the factor right at the end is what i discuss above
        !dEi/dt
        dy(i+n) = (b*C(i)*S(i))-(k+a)*(E(i))-(w)*E(i)+w*(1.d0/dble(dot_product(qnet,A(:,i))))*(dot_product(qnet*A(:,i),E))
        !dCi/dt
        dy(i+(2*n)) = a*E(i)-(g+k)*C(i)-w*C(i)+w*(1.d0/dble(dot_product(qnet,A(:,i))))*(dot_product(qnet*A(:,i),C))
    END DO
    
    

end subroutine rhs



subroutine rhs_omp(n,y,t,a,b0,b1,g,k,w,numthreads,qnet,A,dy)
    implicit none
    !Return RHS of network flu model, parallelized with OpenMP
    !input: 
    !n: total number of nodes
    !y: S,E,C
    !t: time
    !a,b0,b1,g,k,w: model parameters
    !numthreads: the parallel regions should use numthreads threads
    !output: dy, RHS
    integer, intent(in) :: n
    real(kind=8), dimension(n*3),intent(in) :: y
    real(kind=8), intent(in) :: t,a,b0,b1,g,k,w
    real(kind=8), dimension(n*3), intent(out) :: dy
    real(kind=8), dimension(n),intent(in) :: qnet
    real(kind=8), dimension(n,n),intent(in) :: A
    integer :: i
    
    
    

end subroutine rhs_omp



end subroutine rhs_omp
