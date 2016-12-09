from networktest import network
from rwnettest import rwmodule


N0 = 5
L = 2
Nt = 4


qmax,qnet,enet = network.generate(N0,L,Nt)
print enet
alist1, alist2 = network.adjacency_list(qnet,enet)
print alist1
print alist2


Ntime = 5,

Nm = 6
X0 = 2
isample = 1


print rwmodule.rwnet(Ntime,Nm,X0,N0,L,Nt,isample)
