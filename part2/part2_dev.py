"""Project, part 2"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from n1 import network as net
#from p1 import flunet as fn
from time import clock,time


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
    InitialConditions = np.zeros((3*N,1))
    InitialConditions[:N,0] = float(1)
    InitialConditions[InfectedNode,0] = 0.1
    InitialConditions[InfectedNode+N,0] = 0.05
    InitialConditions[InfectedNode+2*N,0] = 0.05
    
    if pflag == True:
        #find the adjacency matrix, A
        A = net.adjacency_matrix(N,enet)
        P = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                P[i,j] = float(qnet[i]*A[i,j])/float(np.dot(qnet,A[:,j]))
                
        return InitialConditions, InfectedNode, P
        
    return InitialConditions, InfectedNode
        
        


    
def solveFluNet(T,Ntime,a,b0,b1,g,k,w,y0,N0,L,Nt):
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
    """
    InitialConditions,InfectedNode,P = initialize(N0,L,Nt,True)
    N = N0+Nt
    
    t=np.linspace(0,T,Ntime)
    
    
    #add input variables to RHS functions if needed
    def RHSnet(y,t,a,b0,b1,g,k,w,N,i):
        """RHS used by odeint to solve Flu model
        requires i, the index of the node being solved
        and N, the total number of nodes
        """
        #Initialise initial conditions
        S = y[0:N]
        E = y[N:2*N]
        C = y[(2*N):3*N]
        

        b = b0 + b1*(1.0+np.cos(2.0*np.pi*t))
        
        dS = k*(1-S[i])-S[i]*((b*C[i])+w)+(w*(np.dot(S,P[:,i])))
        dE = (b*C[i]*S[i]) - E[i]*(k+a+w)+(w*(np.dot(E,P[:,i])))
        dC = (a*E[i])-(C[i]*(g+k+w))+(w*(np.dot(C,P[:,i])))
        
        dy = dS, dE, dC
        return dy
    for i in range(N):
        myInit = [InitialConditions[i],InitialConditions[N+i],InitialConditions[2*N+i]]    
        Y = odeint(myInit,InitialConditions,t,args=(a,b0,b1,g,k,w,N,i))
        
        
    
        
     
    
    def RHSnetF(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model"
        Calculations carried out by fn.rhs
        """
        
        return dy
        
    def RHSnetFomp(y,t,a,b0,b1,g,k,w):
        """RHS used by odeint to solve Flu model
        Calculations carried out by fn.rhs_omp
        """

        return dy

    #Add code here and to RHS functions above to simulate network flu model

    return t,S,E,C


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

    return Cmax,Tmax,Nmax
    

def visualize(enet,C,threshold):
    """Optional, not for assessment: Create figure showing nodes with C > threshold.
    Use crude network visualization algorithm from homework 3
    to display nodes. Contagious nodes should be red circles, all
    others should be black"""
    return None


def performance():
    """function to analyze performance of python, fortran, and fortran+omp approaches
        Add input variables as needed, add comment clearly describing the functionality
        of this function including figures that are generated and trends shown in figures
        that you have submitted
    """


if __name__ == '__main__':            
   a,b0,b1,g,k,w = 45.6,750.0,0.5,73.0,1.0,0.1
   InitialConditions, InfectedNode, P = initialize(5,2,3,True)
   #print(InitialConditions)
   #print(InfectedNode)
   print(P)
   
