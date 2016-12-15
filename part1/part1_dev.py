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
    tan = time.time()
    for i in range(Ntime+1):
        n = np.bincount(X[i,:])[1:]/float(m)
        FsecondPrime[i,:len(n)] = n
    
    print Fprime[:,0]
    print FsecondPrime[:,0]
    tam = time.time()
    
    
    print "rwmodule time is %f" % (t3-t2)
    print "making F time is %f" % (t1-t0)
    print "second making F time is %f" % (tam-tan)
    
    tkeval = time.time()
    if display == True:
        time2 = range(1,Ntime+2)
        maxArray = np.zeros(Ntime+1)
        for i in range(Ntime+1):
            maxArray[i] = int(np.argmax(Fprime[i,:])+1)
        
        plt.figure()
        plt.xlabel('t')
        plt.ylabel('node with greatest number of walkers')
        plt.title('The nodes with the greatest number of walkers at each step')
        plt.xlim([0,Ntime+2])
        plt.scatter(time2,maxArray)
        plt.show()
    tnassim = time.time()
    print "time plotting is %f" % (tnassim-tkeval)
    
    return Fprime
            
    
 


def convergence_rnet(Ntime, m, X0, N0, L, Nt, display):
    """Input variables:
    Ntime: number of time steps
    m: number of walks
    X0: initial node, (node with maximum degree if X0=0)
    N0,L,Nt: recursive network parameters
    """
    
    Fprime = analyze_rnet(Ntime,m,X0,N0,L,Nt,False)
    #F = F[Ntime,X0]
    F = Fprime[Ntime-1,X0-1]
    
    Fprimevals = np.zeros(m)
    for i in range(1,m+1):
        Fprime = analyze_rnet(Ntime,i,X0,N0,L,Nt,False)
        F = Fprime[Ntime-1,X0-1]
        Fprimevals[i-1] = F
    
    plt.figure()
    plt.plot(range(m),Fprimevals)
    plt.show()
    	


if __name__ == '__main__':
    # add code here to call functions and generate figures
    #Ntime,m,X0,N0,L,Nt,display = (100,500,0,5,2,5,True)
    N0,L,Nt = (5,2,200)
    
    #convergence_rnet(Ntime,m,X0,N0,L,Nt,display)
    Ntime,m = 100,10
    X0 = 0
    display = True
#    convergence_rnet(Ntime,m,X0,N0,L,Nt,display)    
    analyze_rnet(Ntime,m,X0,N0,L,Nt,True)
    #second part
    
    
