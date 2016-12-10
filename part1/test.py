from networktest import network
from rwnettest import rwmodule


N0 = 5
L = 2
Nt = 3


qmax,qnet,enet = network.generate(N0,L,Nt)

alist1, alist2 = network.adjacency_list(qnet,enet)




Ntime = 10

Nm = 5
X0 = 2
isample = 1



X,XM = rwmodule.rwnet(Ntime,Nm,X0,N0,L,Nt,isample)


print X
print XM




