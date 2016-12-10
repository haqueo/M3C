from networktest import network
from rwnettest import rwmodule
import time

N0 = 5
L = 2
Nt = 3


qmax,qnet,enet = network.generate(N0,L,Nt)

alist1, alist2 = network.adjacency_list(qnet,enet)




Ntime = 10

Nm = 5
X0 = 0
isample = 3






start1 = time.time()
X,XM = rwmodule.rwnet_omp(Ntime,Nm,X0,N0,L,Nt,isample,2)
start2 = time.time()
print start2-start1 


start3 = time.time()
Y,YM = rwmodule.rwnet(Ntime,Nm,X0,N0,L,Nt,isample)
start4 = time.time()
print start4-start3
