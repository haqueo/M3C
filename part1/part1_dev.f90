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
    integer, dimension(Ntime+1,Nm), intent(out) :: X,XM
    integer :: i1,qmax, node,lowerBound,upperBound,degNode,increment, counter,j
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet
    integer, dimension((N0+L*Nt)*2) :: alist1
    integer, dimension(N0+Nt) :: alist2
    real(kind=8) :: u,current

    
     
    
    !A walker is at node Xi at time, ti
    !at time ti+1 it will move to one of the adjacent nodes
    
    !If X0 = 0, the random walks should begin at the node with the 
    !highest degree
    
    !XM(t) is the fraction of walkers that are at the initial node.
    
    !add code to compute XM every isample time steps
    
    !find the index in qnet with value qmax
    
    
    !generate the network
    call generate(N0,L,Nt,qmax,qnet,enet)
    
    !convert qnet, enet into alist1, alist2
    call adjacency_list(qnet,enet,alist1,alist2)
    
    !Define the initial node
    node = X0
    IF (node .eq. 0) THEN
        do i1 = 1, N0+Nt
            if (qnet(i1) .eq. qmax) then
                node = i1
                exit
            end if
        end do 
    END IF
    
    !complete a single walk
    
    !initialise first row of X
    X(1,:) = node
    
    DO j = 1,Nm
    
    
        DO i1 = 1,Ntime
            !!!GIVEN SUBLIST FUNCTION(node) !!!!
            
            !! GENERATE THE SUBLIST!!!
            !! lower bound is alist1(node) 
            !! upper bound is alist2(node+1)
            !!unless node is the final node. then it's alist1(node) till 
            !the last value of alist2
            
            lowerBound = alist1(node)
            
            
            IF (node .eq. N0+NT) THEN
                    upperBound = size(alist2) ! <- we know what this size is
            ELSE 
                    upperBound = alist1(node+1)
            END IF
            
            
            !now we have our sublist.
            
            degNode = size(alist2(lowerBound:upperBound))
            increment = dble(1)/(dble(degNode))
            
            call random_number(u)
            counter = 0
            current = 0
            
            DO WHILE (current .lt. u)
            
                current = current + increment
                counter = counter + 1
                
            END DO
            
            node = alist2(lowerBound + counter -1)
            !this is the next node
            !so add it to 
            
            X(1+i1,j) = node
        
        
        END DO
        
    END DO
    
    


end subroutine rwnet


subroutine rwnet_omp(Ntime,Nm,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
	implicit none
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample,numthreads
    integer, dimension(Ntime+1,Nm), intent(out) :: X,XM



end subroutine rwnet_omp


end module rwmodule
