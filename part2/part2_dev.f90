module flunet
        use omp_lib
	implicit none
        !add variables as needed
	save
	contains


subroutine rhs(n,y,t,a,b0,b1,g,k,w,P,dy)
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
    real(kind=8), dimension(n,n), intent(in) :: P
    real(kind=8), dimension(3*n), intent(out) :: dy
    integer :: i1
    real(kind=8) :: b
    real(kind=8), dimension(n) :: S,E,C
    real(kind=8), dimension(n) :: dS,dE,dC
    
    !print *, "n is", n
    !print *, "y is", y
    !print *, "t,a,b0,b1 is", t,a,b0,b1
    !print *, "g,k,w is", g,k,w
    !print *, "P is", P
    
    
    
    S = y(1:n) ! S(i) = y(i)
    E = y(n+1:2*n) ! E(i) = y(i+N)
    C = y(2*n+1:3*n) 
    ! C(i) = y(i+(2*N))
    
    b = b0 + b1*(1+cos(dble(8)*atan(dble(1))*t))

    !S[:] = y(1:n)
    !E[:] = y(n+1:2*n)
    !C[:] = y((2*n)+1:3*n)
    
    
    
    DO i1 = 1,N
        !dy(i1) = (k*(1-S(i1))) - (S(i1)*(b*C(i1)))+(w*(dot_product(P(i1,:),S)))-w*S(i1) ! this doesnt work
        !!DONT USE dy(i1+N) = (b*C(i1)*S(i1)) - (E(i1)*(k+a))+(w*(dot_product(E,P(i1,:)))) - w*E(i1)
        !dy(i1+N) = b*C(i1)*S(i1)-(k+a)*E(i1)+w*dot_product(P(i1,:),E)-w*E(i1) !only E works
        !!DONT USE dy(i1+(2*N)) = (a*E(i1))-(C(i1)*(g+k))+(w*(dot_product(C,P(i1,:)))) - w*C(i1)
        !dy(i1+2*N) = a*E(i1) - (g+k)*C(i1) + w*dot_product(P(i1,:),C)-w*C(i1) ! this doesnt work
        
        dS(i1) = (k*(1-S(i1))) - (S(i1)*(b*C(i1)))+(w*(dot_product(P(i1,:),S)))-w*S(i1)
        dE(i1) = b*C(i1)*S(i1)-(k+a)*E(i1)+w*dot_product(P(i1,:),E)-w*E(i1)
        dC(i1) = a*E(i1) - (g+k)*C(i1) + w*dot_product(P(i1,:),C)-w*C(i1)
       
    END DO
    
    dy(1:n) = ds(1:n)
    dy(n+1:2*n) = de(1:n)
    dy(2*n+1:3*n) = dc(1:n)
    !... the column with the initial values is messed up
    
    !2nd, 5th and 8th columns work
    !need use P values without creating an NxN matrix
    
    
    !Choices: (i) add P as an input variable and use it
    !         (ii) add qnet, A as input variables then use them
    ! Considering we need qnet, A to make P anyway - (ii) is probably cheaper
    !but this still uses an NxN matrix..
    !yeah, this is cheaper since we can use vectorisation.
    !no, there's no way to get qnet, and A without making initialize() return it :(
    
    
    
    
    !from python:
    !P[i,j] = float(qnet[i]*A[i,j])/float(np.dot(qnet,A[:,j]))
    !fortran equivalent:
    !P(i,j) = dble(qnet(i)*A(i,j))/float(dot_product(qnet,A(:,j)))
    
    
    !Now, what we want is terms of the form: P(:,i)*dotproduct*S
    !P(:,i)*DOTPRODUCT*S = [qnet(1)*A(1,i)/dot_product(qnet,A(:,i)),...,qnet(N)*A(N,i)/dot_product(qnet,A(:,N))]*dotproduct*S
    !                    = (1/dot_product(qnet,A(:,i))[qnet(1)*A(1,i),qnet(2)*A(1,i),..,qnet(N)*A(N,i)]*dotproduct*S
    !                    = (1/dot_product(qnet,A(:,i))*(qnet*elementwisemultiplication*A(:,i))*dotproduct*S
    !similar for E,C.
    
    !S(i) = y(i)
    !E(i) = y(n+i)
    !C(i) = y((2*n)+i)
    
   ! 
   !this was my first approach, it's fast but it needs qnet and Anet separately
   !which I can't get unless initialize returns them
    !DO i = 1,n
     !   !dSi/dt 
      !  dy(i) = k*(1-y(i))-(b*y((2*n)+i)*y(i))-(w*y(i))+ &
       ! & w*(1.d0/dble(dot_product(qnet,Anet(:,i))))*(dot_product(qnet*Anet(:,i),y(1:n)))
        !!dEi/dt
        !dy(i+n) = (b*y((2*n)+i)*y(i))-(k+a)*(y(n+i))-(w)*y(n+i)+ &
        !& w*(1.d0/dble(dot_product(qnet,Anet(:,i))))*(dot_product(qnet*Anet(:,i),y(n+1:2*n)))
        !dCi/dt
        !dy(i+(2*n)) = a*y(n+i)-(g+k)*y((2*n)+i)-w*y((2*n)+i)+&
        !& w*(1.d0/dble(dot_product(qnet,Anet(:,i))))*(dot_product(qnet*Anet(:,i),y((2*n)+1:3*n)))
!    END DO
    
    
    !!! NEW APPROACH!!! 
    !do exactly the same as in python
    
    




end subroutine rhs



subroutine rhs_omp(n,y,t,a,b0,b1,g,k,w,numthreads,qnet,Anet,dy)
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
    real(kind=8), dimension(n,n),intent(in) :: Anet
    real(kind=8) :: b
    integer :: i,numThreads
    !$ call omp_set_num_threads(numThreads)    
    
    b = b0 + b1*(1.d0+cos(2.0*4.0*atan(1.0)*t))
    !$OMP parallel do
    DO i = 1,n
        !dSi/dt 
        dy(i) = k*(1-y(i))-(b*y((2*n)+i)*y(i))-(w*y(i))+w*(1.d0/dble(dot_product(qnet,Anet(:,i)))) &
        & *(dot_product(qnet*Anet(:,i),y(1:n))) 
        !dEi/dt
        dy(i+n) = (b*y((2*n)+i)*y(i))-(k+a)*(y(n+i))-(w)*y(n+i)+&
        & w*(1.d0/dble(dot_product(qnet,Anet(:,i))))*(dot_product(qnet*Anet(:,i),y(n+1:2*n)))
        !dCi/dt
        dy(i+(2*n)) = a*y(n+i)-(g+k)*y((2*n)+i)-w*y((2*n)+i)+ &
        & w*(1.d0/dble(dot_product(qnet,Anet(:,i))))*(dot_product(qnet*Anet(:,i),y((2*n)+1:3*n)))
    END DO
    !$OMP end parallel do


end subroutine rhs_omp

end module flunet