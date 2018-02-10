This is the code used in Wieland et al., Phys. Rev. E 91, 060801(R) – 2015. It is a Gillespie implementation of the adaptive contact process on complex networks proposed by Gross et al., Phys. Rev. Lett. 96, 208701 – 2006. In this toy model, a dynamic process on a network - disease spreading in a human population - coevolves with that network's topology through a simple rewiring mechanism which mimicks a change in contact patterns during the epidemic. More specifically, the model features three elementary processes: 

(1) disease transmission with rate P along a randomly selected SI-link connecting a susceptible (S-node) and infected individual (I-node): SI->II

(2) recovery with rate R of a randomly selected I-node: I->S

(3) rewiring with rate W of the infected end of a randomly selected SI-link to a randomly selected S-node: SI+S->SS+I

In the process, no multiple or self-connections are allowed. For large chunks of parameter space, this coevolutionary dynamics is observed to yield a simple active phase where infection, recovery and rewiring processes perpertually reshape network topology and node states on a *local* level, while *global* measures of node dynamics and network structure settle down to steady state values. These measures are, among others, the percentage [I] of I-nodes in the network, the per-capita number [SI] of links connecting S-nodes with I-nodes, and the joint degree distributions P{S,I}(x,y), where P{A}(x,y) is the fraction of nodes of state A\in{S,I} with x S-neighbors and y I-neighbors.

Averaging over multiple realizations (~100) of the stochastic process yields decent statistics. In each realization, on starts with a randomly configured Erdos-Renyi graph with a fixed number of initial I-nodes and then lets the coevolutionary dynamics run its course until a steady state is reached (rule of thumb: computed time ~ 1000). The code uses multithreading via OpenMP, maximizing the number of threads. The output file contains, in its first line (commented out), the averaged values of [I] and [SI], as well at the computation time. The subsequent lines sample steady-state P{S,I}(x,y) and have 4 columns: x, y, P{S}(x,y), P{I}(x,y).
