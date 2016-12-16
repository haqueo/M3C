"""Project, part 2"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from n1 import network as net
from p1 import flunet as fn
import time


def initialize(N0,L,Nt,pflag):
    """Generate network, and initial conditions
    If pflag is true, construct and return transport matrix
    """
    #Generate network
    qmax,qnet,enet = net.generate(N0,L,Nt)
    #THERE ARE N = N0+NT NODES IN TOTAL
    N = N0+Nt
    
    #locate the first node with the highest degree
    for i in range(N0+Nt):
        if (qnet[i] == qmax):
            InfectedNode = i #The actual value of InfectedNode is i+1, this is just its index
            break
    
    
    
    
    #Setup the initial conditions
    InitialConditions = np.zeros((3*N,1)) #It is 3Nx1
    InitialConditions[:N,0] = float(1) #first N are all 1, rest are 0
    #apart from the infected node:
    InitialConditions[InfectedNode,0] = float(0.1)
    InitialConditions[InfectedNode+N,0] = 0.05
    InitialConditions[InfectedNode+2*N,0] = 0.05
    
    
    
    if pflag == True:
        #find the adjacency matrix, A
        A = net.adjacency_matrix(N0+Nt,enet)
        
        P = np.zeros((N,N))
        
        scalar = np.dot(qnet,A).astype(float)
        #for i in range(N):
        #    for j in range(N):
        #        P[i,j] = float(qnet[i]*A[i,j])/float(np.dot(qnet,A[:,j]))
        for i in range(N):
            P[i,:] = (float(qnet[i])*A[i,:].astype(float))/scalar
                
        return InitialConditions[:,0], InfectedNode, P
        
    return InitialConditions[:,0], InfectedNode
        


        


    
def solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,N0,L,Nt,P,verFlag):
    """Simulate network flu model at Ntime equispaced times
    between 0 and T (inclusive). Add additional input variables
    as needed clearly commenting what they are and how they are  
    used
    additional input:
        P - the transport matrix assosciated with the generated network,
        used to calculate RHS
        N0 - initial number of nodes 
        L - number of links added at each timestep
        Nt - total number of timesteps
        verFlag - tells you which version of RHS to use
                - "python", "fortran", "omp"
    """
    
    
    
    #add input variables to RHS functions if needed
    def RHSnet(y,t,a,b0,b1,g,k,w,N,P):
        """RHS used by odeint to solve Flu model
        requires i, the index of the node being solved
        and N, the total number of nodes
        """
        
        #y is really the initial condition, this needs to change
        #Split initial conditions into S,E,C
        
        
        
        
        #beta (b) definition
        b = b0 + (b1*(float(1)+np.cos(2.0*(np.pi)*float(t))))
        
    
        #Here I have expanded the sum within dS/dt, dE/dt and dC/dt given,
        # and used the fact that the sum of Pij over i = 1
        #
        #dS = k*(1-S[i])- S[i]*((b*C[i])+w)+(w*(np.dot(S[:,0],P[:,i])))
        #dE = (b*C[i]*S[i]) - E[i]*(k+a+w)+(w*(np.dot(E[:,0],P[:,i])))
        #dC = (a*E[i])-(C[i]*(g+k+w))+(w*(np.dot(C[:,0],P[:,i])))
        #
        #
        dyprime = np.zeros((3*N,1))
        S = np.zeros(N)
        E = np.zeros(N)
        C = np.zeros(N)
        S[:] = y[0:N]
        E[:] = y[N:2*N]
        C[:] = y[2*N:3*N]
        
        
        
        for i in range(N):
            
            #dyprime[i,0] = k*(1-y[i])- y[i]*((b*y[(2*N)+i])+w)+(w*(np.dot(y[0:N],P[i,:])))
            dyprime[i,0] = k*(1-S[i]) - b*C[i]*S[i] + w*np.dot(P[i,:],S) - w*S[i]
            dyprime[i+N,0] = (b*y[(2*N)+i]*y[i]) - y[N+i]*(k+a+w)+(w*(np.dot(y[N:2*N],P[i,:]))) 
            dyprime[i+(2*N),0] = (a*y[N+i])-(y[(2*N)+i]*(g+k+w))+(w*(np.dot(y[2*N:3*N],P[i,:])))
            
        
        
        
        
        return dyprime[:,0]
    #def RHSnet(y,t,a,b0,b1,g,k,w,N,initialConditions,P):
    
    
    
    #
    #for i in range(N):
    #    myInit = InitialConditions[i][0],InitialConditions[N+i][0],InitialConditions[2*N+i][0]
    #    #for each node, extract its initial condition (could use if statement dependent on if infectedNode or not, but still O(1)
    #    
    #    Y = odeint(RHSnet,myInit,t,args=(a,b0,b1,g,k,w,N,i,InitialConditions)) #i is the node under consideration
    #    
    #    S_sol[:,i] = Y[:,0]
    #    E_sol[:,i] = Y[:,1]
    #    C_sol[:,i] = Y[:,2]
        
    
    #solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,N0,L,Nt):
        
        
    
        
     
    
    def RHSnetF(y,t,a,b0,b1,g,k,w,N,P):
        """RHS used by odeint to solve Flu model"
        Calculations carried out by fn.rhs
        """
        
       
        #!rhs(n,y,t,a,b0,b1,g,k,w,dy,qnet,Anet)
        dy = fn.rhs(y,t,a,b0,b1,g,k,w,P,N)
        #for some reason, this takes N as the last value?!
        
        return dy
        
        
    
    
    
    def RHSnetFomp(y,t,a,b0,b1,g,k,w,N,P):
        """RHS used by odeint to solve Flu model
        Calculations carried out by fn.rhs_omp
        """
        dy = fn.rhs_omp(y,t,a,b0,b1,g,k,w,P,2,N)
        return dy
        
        
    #there are N0+Nt nodes in total
    N = N0+Nt
    
    #Get the initial conditions
    
    
    #Initialise the time
    tprime = np.linspace(float(0),T,Ntime)
   
   
    
    #Add code here and to RHS functions above to simulate network flu model

    if (verFlag == 'python'):
        Y = odeint(RHSnet,y0,tprime,args=(a,b0,b1,g,k,w,N,P)) 
        return tprime, Y[:,0:N],Y[:,N:2*N],Y[:,(2*N):3*N]
      
    
    elif (verFlag == 'fortran'):
        Y = odeint(RHSnetF,y0,tprime,args=(a,b0,b1,g,k,w,N,P)) 
        return tprime, Y[:,0:N],Y[:,N:2*N],Y[:,(2*N):3*N]
       
    
    elif (verFlag == 'omp'):
        Y = odeint(RHSnetFomp,y0,tprime,args=(a,b0,b1,g,k,w,N,P)) 
        return tprime, Y[:,0:N],Y[:,N:2*N],Y[:,(2*N):3*N]
        

def analyze(N0,L,Nt,T,Ntime,a,b0,b1,g,k,threshold,warray,display=False):
    """analyze influence of omega on: 
    1. maximum contagious fraction across all nodes, Cmax
    2. time to maximum, Tmax
    3. maximum fraction of nodes with C > threshold, Nmax    
    Input: N0,L,Nt: recursive network parameters 
           T,Ntime,a,b0,b1,g,k: input for solveFluNet
           threshold: use to compute  Nmax
           warray: contains values of omega
    """
    
    Wlen = len(warray)
    Cmaxarray = np.zeros((Wlen,1))
    Tmaxarray = np.zeros((Wlen,1))
    Nmaxarray = np.zeros((Wlen,1))
    InitialConditions, InfectedNode, P = initialize(N0,L,Nt,True)
    N = N0+Nt
    
    
    for i in range(Wlen):
        
        #FOR EVERY W, SOLVE FLU NET 
        t,S,E,C = solveFluNet(T,Ntime,a,b0,b1,g,k,warray[i],InitialConditions,N0,L,Nt,P,'omp')
        
        Tmaxarray[i,0] = t[np.argmax(C,0)[0]]
        
        arrayNs = np.zeros((Ntime,1))
        
        for j in range(Ntime):
            
            #TAKE VALUES OF C OVER THRESHOLD
            arrayNs[j,0] = (C[j,:] > threshold).sum()
            
            #FIND MAX X
            Cmaxarray[i] = np.amax(C)
            
            #NUMBER OF Cs over THRESHOLD
            Nmaxarray[i] = max(arrayNs)/float(N)
            
    if display == True:
        
        plt.figure()
        plt.title('Omar Haque, analyze. Graph of Nmax against w')
        plt.xlabel('w')
        plt.ylabel('Nmax')
        plt.plot(warray,Nmaxarray)
        plt.show()
        
        plt.figure()
        plt.title('Omar Haque, analyze. Graph of Cmax against w')
        plt.xlabel('w')
        plt.ylabel('Cmax')
        plt.plot(warray,Cmaxarray)
        plt.show()
        
        
        plt.figure()
        plt.title('Omar Haque, analyze. Graph of Tmax against w')
        plt.xlabel('w')
        plt.ylabel('Tmax')
        plt.plot(warray,Tmaxarray)
        plt.show()
        

    
    

    return Cmaxarray,Tmaxarray,Nmaxarray
    

def visualize(enet,C,threshold):
    """Optional, not for assessment: Create figure showing nodes with C > threshold.
    Use crude network visualization algorithm from homework 3
    to display nodes. Contagious nodes should be red circles, all
    others should be black"""
    return None


def performance(N0,L,Nt,T,Ntime,a,b0,b1,g,k,w):
    """function to analyze performance of python, fortran, and fortran+omp approaches
        Add input variables as needed, add comment clearly describing the functionality
        of this function including figures that are generated and trends shown in figures
        that you have submitted
        
        This function plots the time taken to complete solveFluNet, against N (total number of nodes), for the three RHS functions.
        The first (performance1) shows how the two Fortran functions are far superior than the python version, and it seems to grow linearly with the problem size.
        
        The second (performance2) gives a log-log plot of the Fortran RHS with the Fortran+OpenMP RHS functions, against N. It is clear that there does not seem to be much of a difference between them, although
            the openMP version does seem to take over for larger N
        
        
        
    """
    N = N0+Nt
    
    Nvalues = np.arange((N))
    pythonTimes = np.zeros(N)
    fortranTimes = np.zeros(N)
    openMPTimes = np.zeros(N)
    
    for i in Nvalues:
        
        #Generate initial values
        
        InitialConditions, InfectedNode, P = initialize(N0,L,Nvalues[i],True)    
        pythonTime1 = time.time()
        #solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,N0,L,Nt,P,verFlag):
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,InitialConditions,N0,L,Nvalues[i],P,'python')
        pythonTime2 = time.time()
        
        fortranTime1 = time.time()
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,InitialConditions,N0,L,Nvalues[i],P,'fortran')
        fortranTime2 = time.time()
        
        openMPtime1 = time.time()
        solveFluNet(T,Ntime,a,b0,b1,g,k,w,InitialConditions,N0,L,Nvalues[i],P,'omp')
        openMPtime2 = time.time()
        #solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,N0,L,Nt,verFlag):
        
        pythonTimes[i] = pythonTime2-pythonTime1
        fortranTimes[i] = fortranTime2-fortranTime1
        openMPTimes[i] = openMPtime2-openMPtime1
        
        
    plt.figure()
    plt.title('Omar Haque, performance. Comparison of different RHS methods')
    plt.plot(Nvalues,pythonTimes,label='RHS from python')
    plt.plot(Nvalues,fortranTimes,label='RHS from fortran')
    plt.plot(Nvalues, openMPTimes,label='RHS using openMP')
    plt.xlabel('N')
    plt.ylabel('time taken for solveFluNet()')
    plt.legend(loc = 'best')
    plt.show()
    
    plt.figure()
    plt.title('Omar Haque, performance. Comparison of Fortran vs Fortran+OpenMP (logplot)')
    plt.plot(np.log(Nvalues),np.log(fortranTimes),label='RHS from fortran')
    plt.plot(np.log(Nvalues), np.log(openMPTimes),label='RHS using openMP')
    plt.xlabel('log(N)')
    plt.ylabel('log(time taken for solveFluNet)')
    plt.legend(loc = 'best')
    plt.show()

        
    
    


if __name__ == '__main__':            
   a,b0,b1,g,k,w = 45.6,750.0,0.5,73.0,1.0,0.1
   #CODE FOR P21,P22,P23
   N0,L,Nt = 5,2,500
   T = 2
   threshold = 0.1
   warray = [0.0,0.01,0.1,0.2,0.5,1.0]
   display = True
   #def analyze(N0,L,Nt,T,Ntime,a,b0,b1,g,k,threshold,warray,display=False):
   analyze(N0,L,Nt,T,1000,a,b0,b1,1000,k,threshold,warray,True)

   
   #CODE FOR PERFORMANCE
   #I have chosen Nt = 200 as it covers the section where Fortran+OpenMP is faster than Fortran.
   #InitialConditions, InfectedNode, P = initialize(5,3,3,True)
   #N0,L,Nt = 5,2,200
   #T = 2
   #threshold = 0.1
   #Ntime = 1000
   #warray = [0,np.power(10,-2),0.1,0.2,0.5,1.0]
   #performance(N0,L,Nt,T,Ntime,a,b0,b1,g,k,w)
   

   
