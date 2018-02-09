#!/bin/bash

g++ -Wall -g -std=c++11 -O3 -o adaptiveContact adaptiveContact.cpp -lm -lgsl -lgslcblas -fopenmp
./sis 1000 5000 500 1 1 1 1000 100 1000 50 'test.dat'

#./sis N K I P R W TMAX RUNS EMAX KMAX PATH
#
#----model parameters--------------------
#N=number of network nodes
#K=number of network links
#I=inital number of I-nodes
#P=infection rate
#R=recovery rate
#W=rewiring rate
#
#----simulation parameters--------------------
#TMAX=simulated time
#RUNS=number of runs
#EMAX=number of maxixum tries for rewiring procedure
#KMAX=maximum considered degree
#PATH=path for data files

