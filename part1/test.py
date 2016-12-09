from networktest import network

qmax,qnet,enet = network.generate(12,4,8)
print enet
alist1, alist2 = network.adjacency_list(qnet,enet)
print alist1

