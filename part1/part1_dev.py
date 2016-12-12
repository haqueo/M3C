
import numpy as np
import matplotlib.pyplot as plt
import rw #Assumes rwmodule has been compiled with f2py to produce rw.so


def analyze_rnet(Ntime,m,X0,N0,L,Nt,display):
	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    #The function, analyze_rnet should call rwnet or rwnet_omp and analyze the results. The routine should compute and return, F(t,n), the fraction of the m walkers at node n at time, t. When display is true, a figure should be created which plots the node with the greatest number of walkers at each step.







def convergence_rnet(Ntime,m,X0,N0,L,Nt,display):
    	"""Input variables:
	Ntime: number of time steps
    	m: number of walks
    	X0: initial node, (node with maximum degree if X0=0)
    	N0,L,Nt: recursive network parameters
    	"""
    	
    	
    	
if __name__== '__main__':
#add code here to call functions and generate figures

