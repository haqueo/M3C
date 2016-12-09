!Project part 1
module rwmodule
    use network

contains



subroutine rwnet(Ntime,m,X0,N0,L,Nt,isample,X,XM)
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
    integer :: i1,qmax
    integer, dimension(N0+Nt) :: qnet
    integer, dimension(N0+L*Nt,2) :: enet 
    
    !A walker is at node Xi at time, ti
    !at time ti+1 it will move to one of the adjacent nodes
    
    !If X0 = 0, the random walks should begin at the node with the 
    !highest degree
    
    !XM(t) is the fraction of walkers that are at the initial node.
    
    !add code to compute XM every isample time steps
    
    call generate(N0,L,Nt,qmax,qnet,enet)
    


end subroutine rwnet


subroutine rwnet_omp(Ntime,m,X0,N0,L,Nt,isample,numthreads,X,XM)
    !parallelized version of rwnet, parallel regions should
    !use numthreads threads
	implicit none
    integer, intent(in) :: Ntime,Nm,N0,L,Nt,X0,isample,numthreads
    integer, dimension(Ntime+1,Nm), intent(out) :: X



end subroutine rwnet


end module rwmodule
