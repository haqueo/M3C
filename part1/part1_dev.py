import numpy as np
import matplotlib.pyplot as plt
from rwnettest import rwmodule 






# Assumes rwmodule has been compiled with f2py to produce rw.so
def counter(number,row):
    """Counts the number of occurences of an element
        in a given list
    """
    numberCount = 0
    for _,item in enumerate(row):
        if item == number:
            numberCount = numberCount + 1
    return numberCount



def analyze_rnet(Ntime, m, X0, N0, L, Nt, display):
    """Input variables:
	Ntime: number of time steps
    m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	Output: Fprime - an (Ntime+1)x(N0+Nt) matrix
    	where Fprime[t-1,n-1] = F(t,n) from the question
    	"""
    X,XM = rwmodule.rwnet(Ntime,m,X0,N0,L,Nt,Ntime)
    Fprime = np.zeros((Ntime+1,N0+Nt))
    
    for i in range(Ntime+1):
        for j in range(N0+Nt):
            Fprime[i,j] = \
            float(counter(j+1,X[i,:])) \
            /float(m)
            
    print Fprime
            
            
    
            
    print X
 
    
    
    
analyze_rnet(5,6,2,8,2,2,True)

# The function, analyze_rnet should call rwnet or rwnet_omp and analyze the results.
# The routine should compute and return, F(t,n), the fraction of the m walkers at node n at time, t. 
# When display is true, a figure should be created which plots the node with the greatest number of walkers at each step.







def convergence_rnet(Ntime, m, X0, N0, L, Nt, display):
    """Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""


if __name__ == '__main__':
    # add code here to call functions and generate figures
    pass
