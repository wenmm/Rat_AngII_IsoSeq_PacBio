import networkx as nx
import sys

a = []
with open(sys.argv[1],'r') as f:
    G = nx.Graph()
    for line in f:
        table = line.rstrip().split()
        i = table[0]+'|'+table[1]                                  
	j = table[0]+'|'+table[2]        
	G.add_edge(i,j)
    h = [g.nodes() for g in nx.connected_component_subgraphs(G)]
    for i in h:
        print(i)
