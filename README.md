# adaptive-contact-process
stochastic simulation of disease spreading in a population with a simple contact avoidance strategy

This is the code used in Wieland et al.,Phys. Rev. E 91, 060801(R) – 2015. It is a Gillespie implementation of the adaptive contact process on complex networks proposed by Gross et al, Phys. Rev. Lett. 96, 208701 – 2006. In this toy model, a dynamic process on a network - disease spreading in a human population - coevolves with that network's topology through a simple rewiring mechanism which mimicks a change in contact patterns during the epidemic. For large chunks of parameter space, this process is observed to yield a simple active phase where infection, recovery and rewiring processes perpertually reshape network topology and node states on a *local* level, while *global* measures of node dynamics and network structure settle down to steady state values. These measures are, among others, the percentage [I] of infected (I-)nodes in the network, the per-capita number [SI] of links connecting susceptible (S-)nodes with I-nodes, and the joint degree distributions P{S,I}(x,y), where P{A}(x,y) is the fraction of A-nodes with x S-neighbors and y I-neighbors.

The output file contains, in its first line, the averaged values of the stationary densities of infected (I-)nodes and of links connecting susceptible (S-)nodes with I-nodes, as well at the computation time. The subsequent lines sample steady-state P{S,I}(x,y) and have 4 columns: 

1st: x
2st: y
3rd: P{S}(x,y)
4th: P{I}(x,y).

For better averaging, multiple runs are possible. In each run, on starts with a randomly configured Erdos-renyi graph with a given number of initial I-nodes and then lets the coevolutionary dynamics run its course until a steady state is reached (rule of thumb: computed time >= 1000). The code uses multithreading via OpenMP, maximizing the number of threads.  
