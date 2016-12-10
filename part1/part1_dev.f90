!Project part 1
module rwmodule
    use network

contains



subroutine rwnet(Ntime,Nm,X0,N0,L,Nt,isample,X,XM)
    !random walks on a recursive network
    !Input variables:
    !Ntime: number of time steps
    !m: number of walks
    !X0: initial node, (node with maximum degree if X0=0)
    !N0,L,Nt: recursive network parameters
    !Output: X, m random walks on network each with Ntime steps from initial node, X0
    !XM: fraction of walks at initial node, computed every isample timesteps
    implicit none
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample
    integer, dimension(Ntime+1,Nm), intent(out) :: X
    real(kind=8), dimension((Ntime+1)/isample),intent(out) :: XM !<- check this is integer division
    real(kind=8), dimension((Ntime+1)/isample) :: XMtemp !<- check this is integer division
    integer :: i1,qmax, node,lowerBound,upperBound,degNode, counter,j,initialNode,initialCounter
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension(size(enet,1)*2) :: alist1
    integer, dimension(size(qnet)) :: alist2
    real(kind=8) :: u,current, increment

     
    
    
    !generate the network
    call generate(N0,L,Nt,qmax,qnet,enet)
    
    !convert qnet, enet into alist1, alist2
    call adjacency_list(qnet,enet,alist1,alist2)
    
    
    
    !initialise XM and XMtemp
    XM(:) = dble(0)
    XMtemp(:) = dble(0)
    
    !Find the initial node. 
    node = X0
    !If the input is not zero then we just keep X0
        IF (node .eq. 0) THEN
    !If it is zero then we need to find the node with highest degree
    !We could use the inbuilt MAXLOC, but we already have qmax
    !So there's no need to go further than an element with value qmax
            DO i1 = 1,N0+Nt
                IF (qnet(i1) .eq. qmax) THEN
                    node = i1
                    EXIT
                end IF
            END DO
        END IF
    !we save this initial node for later.
    initialNode = node
    
    
    DO j = 1,Nm
    
        !"for every column"
        !set the first value as initial node
        X(1,j) = initialNode
        
        
        !we only need this to contribute to XM if isample is 1
        IF (isample .eq. 1) THEN
            XMtemp(1) = 1
        END IF
        
        node = initialNode
        
        !iterate through the rest of the column
      
        DO i1 = 1,Ntime
            
            
            !this is the lowerbound of possible values of X(i1+1)
            lowerBound = alist2(node)
            
            
            !the upper bound is given by alist2(node+1)
            !UNLESS, we're at the final node.
            
            IF (node .eq. N0+NT) THEN
                    upperBound = size(alist1)+1 ! <- we know what this size is
            ELSE 
                    upperBound = alist2(node+1)
            END IF
            
            
            
            !now we have our sublist.
            
            !IDEA:
            !we have a sublist of possible nodes we can go to
            !sublist = [1,2,3,2,4], say.
            !to randomise, we want the probabilities of hitting any number to
            !be equal. So the Probability of hitting any of these numbers should
            !be 1/5 (1/length(sublist))
            !to generate from this list, we can consider having intervals of length
            !0.2, corresponding to each number. The probability of a random uniform
            !being in any of these intervals is therefore 1/5.
            !but rather than computing all of the intervals, we just increment
            !until we find the interval the uniform has landed in.
            
            
            !initialise the degree of the Node and the increment (probability 
            !of each interval)
            degNode = upperBound-lowerBound
            increment = dble(1)/(dble(degNode))
            
            !generate the random uniform u
            call random_number(u)
            
            !initialise counter and current. counter tells you which interval 
            !you're in, and current tells you the current cumulative probability
            
            counter = 0
            current = 0
            
            DO WHILE (current .lt. u)
            
                current = current + increment
                counter = counter + 1
                
            END DO
            
            !we're now at a place where current is greater than or equal to u
            !so we're in the right interval. We can now add this node to X
            
            node = alist1(lowerBound + counter -1)
            
            
            !ADD THE NODE TO X
            X(1+i1,j) = node
            
           
            
            !if this is a multiple of isample, we need to check if it's an 
            !initial node and add it to XM
            
            IF (int(MOD(i1+1,isample)) .eq. 0) THEN
              
                    IF (node .eq. initialNode) THEN
                     
                        XMtemp((i1+1)/isample) = dble(1)
                    END IF
            END IF
        
        
        END DO
        
        XM = XM+XMtemp
        XMtemp(:) = dble(0)
        
        
    END DO
    
    XM = XM/dble(Nm)
    
    
    
end subroutine rwnet


subroutine rwnet_omp(Ntime,Nm,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
    use omp_lib
    implicit none
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample,numthreads
    integer, dimension(Ntime+1,Nm), intent(out) :: X
    real(kind=8), dimension((Ntime+1)/isample),intent(out) :: XM !<- check this is integer division
    real(kind=8), dimension((Ntime+1)/isample) :: XMtemp !<- check this is integer division
    integer :: i1,qmax, node,lowerBound,upperBound,degNode, counter,j,initialNode,initialCounter
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension(size(enet,1)*2) :: alist1
    integer, dimension(size(qnet)) :: alist2
    real(kind=8) :: u,current, increment


    !$ call omp_set_num_threads(numThreads)
    
    !generate the network
    call generate(N0,L,Nt,qmax,qnet,enet)
    
    !convert qnet, enet into alist1, alist2
    call adjacency_list(qnet,enet,alist1,alist2)
    
    
    
    !initialise XM and XMtemp
    XM(:) = dble(0)
    XMtemp(:) = dble(0)
    
    !Find the initial node. 
    node = X0
    !If the input is not zero then we just keep X0
        IF (node .eq. 0) THEN
    !If it is zero then we need to find the node with highest degree
    !We could use the inbuilt MAXLOC, but we already have qmax
    !So there's no need to go further than an element with value qmax
    
    !A DO loop with an exit, is it worth parallelising?
            DO i1 = 1,N0+Nt
                IF (qnet(i1) .eq. qmax) THEN
                    node = i1
                    EXIT
                end IF
            END DO
        END IF
    !we save this initial node for later.
    initialNode = node
    
    !$OMP parallel do private(XMtemp,lowerBound,upperBound,node,degNode,increment,counter,current,u,i1)
    DO j = 1,Nm
    
        !"for every column"
        !set the first value as initial node
        X(1,j) = initialNode
        
        
        !we only need this to contribute to XM if isample is 1
        IF (isample .eq. 1) THEN
            XMtemp(1) = 1
        END IF
        
        node = initialNode
        
        !iterate through the rest of the column
      
        DO i1 = 1,Ntime
            
            
            !this is the lowerbound of possible values of X(i1+1)
            lowerBound = alist2(node)
            
            
            !the upper bound is given by alist2(node+1)
            !UNLESS, we're at the final node.
            
            IF (node .eq. N0+NT) THEN
                    upperBound = size(alist1)+1 ! <- we know what this size is
            ELSE 
                    upperBound = alist2(node+1)
            END IF
            
            
            
            !now we have our sublist.
            
            !IDEA:
            !we have a sublist of possible nodes we can go to
            !sublist = [1,2,3,2,4], say.
            !to randomise, we want the probabilities of hitting any number to
            !be equal. So the Probability of hitting any of these numbers should
            !be 1/5 (1/length(sublist))
            !to generate from this list, we can consider having intervals of length
            !0.2, corresponding to each number. The probability of a random uniform
            !being in any of these intervals is therefore 1/5.
            !but rather than computing all of the intervals, we just increment
            !until we find the interval the uniform has landed in.
            
            
            !initialise the degree of the Node and the increment (probability 
            !of each interval)
            degNode = upperBound-lowerBound
            increment = dble(1)/(dble(degNode))
            
            !generate the random uniform u
            call random_number(u)
            
            !initialise counter and current. counter tells you which interval 
            !you're in, and current tells you the current cumulative probability
            
            counter = 0
            current = 0
            
            DO WHILE (current .lt. u)
            
                current = current + increment
                counter = counter + 1
                
            END DO
            
            !we're now at a place where current is greater than or equal to u
            !so we're in the right interval. We can now add this node to X
            
            node = alist1(lowerBound + counter -1)
            
            
            !ADD THE NODE TO X
            X(1+i1,j) = node
            
           
            
            !if this is a multiple of isample, we need to check if it's an 
            !initial node and add it to XM
            
            IF (int(MOD(i1+1,isample)) .eq. 0) THEN
              
                    IF (node .eq. initialNode) THEN
                     
                        XMtemp((i1+1)/isample) = dble(1)
                    END IF
            END IF
        
        
        END DO
        
        XM = XM+XMtemp
        XMtemp(:) = dble(0)
        
        
    END DO
    !$OMP end parallel do
    
    XM = XM/dble(Nm)



end subroutine rwnet_omp


end module rwmodule
