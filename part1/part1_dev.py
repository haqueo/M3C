import numpy as np
import matplotlib.pyplot as plt
from networktest import network
from rwnettest import rwmodule 
import time 







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
    t2 = time.time()
    Y = rwmodule.rwnet(Ntime,m,X0,N0,L,Nt,Ntime)
    t3 = time.time()
    
    X = Y[0]
    #Fprime = np.zeros((Ntime+1,N0+Nt))
    #t0 = time.time()
    #for i in range(Ntime+1):
    #    for j in range(N0+Nt):
    #        Fprime[i,j] = \
    #        float(counter(j+1,X[i,:])) \
    #        /float(m)
    #t1 = time.time()
    
    #SECOND TRY AT F
    
    FsecondPrime = np.zeros((Ntime+1,N0+Nt))
   
    for i in range(Ntime+1):
        n = np.bincount(X[i,:])[1:]/float(m)
        FsecondPrime[i,:len(n)] = n
    
    #print Fprime[:,0]
    
    
    if display == True:
        time2 = range(1,Ntime+2)
        maxArray = np.zeros(Ntime+1)
        for i in range(Ntime+1):
            maxArray[i] = int(np.argmax(FsecondPrime[i,:])+1)
        
        plt.figure()
        plt.xlabel('t')
        plt.ylabel('node with greatest number of walkers')
        plt.title('The nodes with the greatest number of walkers at each step')
        plt.xlim([0,Ntime+2])
        plt.scatter(time2,maxArray)
        plt.show()
   
    return FsecondPrime
            
    
 


def convergence_rnet(Ntime, m, X0, N0, L, Nt, display):
    """Input variables:
    Ntime: number of time steps
    m: number of walks
    X0: initial node, (node with maximum degree if X0=0)
    N0,L,Nt: recursive network parameters
    
    Plots a graph of m against the assosciated F(Ntime,X0)
    For large m, we can see how F(Ntime,X0) tails out and converges
    to a single value
    (but it takes a long time to run)
    
    """
    
    
    #Draws a graph of m against F(Ntime,X0)
    Fprime = analyze_rnet(Ntime,m,X0,N0,L,Nt,False)
    #F = F[Ntime,X0]
    F = Fprime[Ntime-1,X0-1]
    
    Fprimevals = np.zeros(m)
    for i in range(1,m+1):
        if (i % 10) == 0:
            
            Fprime = analyze_rnet(Ntime,i,X0,N0,L,Nt,False)
            F = Fprime[Ntime-1,X0-1]
            if F == float(0):
                Fprimevals[i-1] = float('nan')
            else:
                Fprimevals[i-1] = F
    
    
    
    plt.figure()
    plt.title('Variation of F(Ntime,X0) on m')
    plt.xlabel('m')
    
    plt.scatter(range(m),Fprimevals)
    plt.ylabel('F(Ntime,X0)')
    plt.show()
    	


if __name__ == '__main__':
    # add code here to call functions and generate figures
    
    #FIRST GRAPH
    #N0,L,Nt = (5,2,200)
    #Ntime,m = 100000,1000
    #X0 = 0
    #display = True
    #analyze_rnet(Ntime,m,X0,N0,L,Nt,display)
    
    #SECOND GRAPH
    N0,L,Nt = (5,2,200)
    Ntime,m = 100000,500
    X0 = 0
    display = True
    convergence_rnet(Ntime,m,X0,N0,L,Nt,display)
    
    
    
    
    
