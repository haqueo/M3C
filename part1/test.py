from networktest import network
from rwnettest import rwmodule



qmax,qnet,enet = network.generate(5,2,4)
print enet
alist1, alist2 = network.adjacency_list(qnet,enet)
print alist1
print alist2

print len(alist1)