#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include "omp.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class nnClass //object representing next neighbor of node, stored in linked lists in class netNode
{
    public:
        long int nodeID;				//next-neighbor ID
        nnClass(long int nodeIDA)			//constructor
        {
            nodeID=nodeIDA;
        };
        ~nnClass(){};					//destructor
};

class netNode					//object defining a single network node
{
    public:
        int state;				//node state (0=S,1=I)
	std::vector< long int> nnState;		//neighborhood composition (nnState[0/1]: # of S/of I-neighbors)
	std::vector<nnClass> nn;		//linked list of next neighbors
        netNode()				//constructor
	{
		nnState.resize(2);		//two node states for SIS dynamics: S=0, I=1
		nnState[0]=0;nnState[1]=0;	//zero S- and I-neighbors before network construction
		state=0;			//initial node state always S
	};
        ~netNode(){};				//destructor
};

void stateInitialization(std::vector<netNode> &net,long int I0, gsl_rng **rng)//randomly infect IO S-nodes in network
{	
	long int n=net.size();				//total network size, n=N
	std::vector< long int> sPool(n);   		//pool of n=N S-nodes to choose from from                                         

	for(long int j=0;j<n;++j)			//tag available S-nodes; initially, S-Node j is assigned slot j
	        sPool[j]=j;
  
	for(long int iNode=0;iNode<I0;++iNode)		//browsing pool of S-nodes I0 times
    	{
        	long int pickedSlot=gsl_rng_uniform_int(*rng,n-iNode);//randomly pick slot among n-I0 remaining slots
 	        net[sPool[pickedSlot]].state=1;	//turn S-node at chosen slot into I-node
        	sPool[pickedSlot]=sPool[n-1-iNode];//move S-node at last slot of pool to slot of picked S-node
    	};
}

void addLink(std::vector<netNode> &net, long int node1,  long int node2)//adding link between nodes with IDs "node1" and "node2"
{
	net[node1].nn.push_back(nnClass(node2));	//updating linked list of next neighbors of node1
	net[node1].nnState[net[node2].state]++;		//updating neighborhood composition of node1
	net[node2].nn.push_back(nnClass(node1));	//updating linked list of next neighbors of node2
	net[node2].nnState[net[node1].state]++;		//updating neighborhood composition of node2
};

void removeLink(std::vector<netNode> &net, long int node1,  long int node1Index)//removing link between nodes with IDs "node1" and its next neighbor at position "node1Index" on linked list nn
{
	long int node2=net[node1].nn[node1Index].nodeID;		//node2-node1 link shall be removed

//----remove next neighbor "node2" from nn of node1----
	net[node1].nn[node1Index]=std::move(net[node1].nn.back());	//moving node2 to end of nn
	net[node1].nn.pop_back();					//deleting last entry of nn
	net[node1].nnState[net[node2].state]--;				//updating neighborhood composition of node1

//----remove next neighbor "node1" from nn of node2----
	long int node2Index=0;						//node1 at slot "node2Index" of nn list of node2 
	while(net[node2].nn[node2Index].nodeID!=node1)			//determine node2Index
		node2Index++;
	net[node2].nn[node2Index]=std::move(net[node2].nn.back());	//moving node1 to end of nn
	net[node2].nn.pop_back();
	net[node2].nnState[net[node1].state]--;				//updating neighborhood composition of node2		
		
};

void ErdosRenyi(std::vector<netNode> &net,long int k,gsl_rng **rng)//construct random graph with k edges
{
    const  long int n=net.size();			//number of network nodes
    for(long int edge=0;edge<k;++edge)			//k iterations to find two yet unconnected nodes
    {
        bool edgeAlreadyExists=false;			//declaring boolean variable 
        long int node1=-1,node2=-1;			//initializing IDs of node pair
        do
        {
		node1=gsl_rng_uniform_int(*rng,n);	//randomly pick first node "node1"
		node2=gsl_rng_uniform_int(*rng,n);	//randomly pick second node "node2"

		edgeAlreadyExists=false;		//assume node pair is not yet connected
	        if(node1!=node2)			//if node1 and node2 are not identical
        	{
			for (auto &i : net[node1].nn)	//see if node1 and node2 are already connected
			{
        	        	if(i.nodeID==node2)	///if node2 is on nn-list of node 
        	        	{
        	        	        edgeAlreadyExists=true;	//...there is already a link between node1 and node2 
					break;           	//...and we can call off the search for the link
        	        	};
        	        };
        	};
        }
        while(node1==node2 || edgeAlreadyExists==true);	//repeat if selected pair is identical or already connected 
	addLink(net,node1,node2);			//add link between node1 and node2 with time tag t=0
    };
}

void countMoments(std::vector<netNode> &net, long int &varI,  long int &varSI)//counting number of I-nodes and SI-links
{
	varI=0;varSI=0;				//propensity counters
	for(auto & j: net)			//browse network
	{
		if(j.state==0)			//if S-node ...
			varSI+=j.nnState[1];	//...increase count of SI-nodes by number I-neighbors
		else			
		{	
			varSI+=j.nnState[0];	//...increase count of SI-nodes by number S-neighbors
			varI++;			//...increase count of I-nodes
		};
	};
	varSI=varSI/2;				//in loop, we counted SI-links twice - once from each end
} 

void printNetwork(std::vector<netNode> &net)	//visualization tool for troubleshooting with small networks
{
	long int n=net.size();			//network size
	std::cout<<"\n";
	for( long int j=0;j<n;++j)		//loop through network
	{
	        if(net[j].state==1)		
		        std::cout<<j<<"+ ";	//mark I-node with "+"
	        else
		        std::cout<<j<<"- ";	//mark S-node with "-"

		
		for(auto &l : net[j].nn) 			//browse neighborhood of node j
	        {
		        if(net[l.nodeID].state==1)		
        		        std::cout<<l.nodeID<<"+";	//mark neighboring I-node with "+"
        		else
        		        std::cout<<l.nodeID<<"-";	//mark neighboring S-node with "-"
        	};
        	std::cout<<"\n";
	};
}

void recovery(std::vector<netNode> &net,long int &varI,long int &varSI,gsl_rng **rng)//recovery process I->S
{	
    	long int pick=gsl_rng_uniform_int(*rng,varI)+1;	//randomly pick I-node
	long int iCounter=0;				//goal: iCounter=pick 
    	long int candidateNode=-1;			//search variable for I-node

	while(iCounter<pick)				//browse to picked I-node
    	{
		candidateNode++;
		if(net[candidateNode].state==1)
			iCounter++;
	};

	net[candidateNode].state=0;		//turn I-node to S-node	
	varI--;					//update I-node count
	varSI+=net[candidateNode].nnState[1]-net[candidateNode].nnState[0];	//update SI-link count
  	for(auto &i : net[candidateNode].nn)	//update neighborhood composition of neighbors of recovered node:
	{
        	net[i.nodeID].nnState[1]--;	//...one less I-neighbor 
        	net[i.nodeID].nnState[0]++;	//...one more S-neighbor
	};
};

void infection(std::vector<netNode> &net, long int &varI, long int &varSI,gsl_rng **rng)//infection process SI->II
{	
	long int candidateLink=gsl_rng_uniform_int(*rng,2*varSI)+1;	//randomly pick SI-link (links are counted twice)
	long int siCounter=0;						//goal: siCounter=candidateLink
	long int candidateNode=-1;					//search variable for one end node of SI-link
	while(siCounter<candidateLink)					//browse SI-links nodewise
	{
		candidateNode++;					//look for one end of SI-link
		if(net[candidateNode].state==0)				//increase SI-link count 
			siCounter+=net[candidateNode].nnState[1];
		else
			if(net[candidateNode].state==1)		
				siCounter+=net[candidateNode].nnState[0];
	};
	int nodeState=net[candidateNode].state;				//determine state of one SI-link end
	long int candidateNodeIndex=net[candidateNode].nn.size()-1;	//starting index of search for other end of SI-link
	while(siCounter>candidateLink || (nodeState+net[net[candidateNode].nn[candidateNodeIndex].nodeID].state!=1))//browse SI-links linkwise and backwards, as siCounter always overshoots in previous loop
	{
		if(nodeState+net[net[candidateNode].nn[candidateNodeIndex].nodeID].state==1)
			siCounter--;
		candidateNodeIndex--;
	};

	long int sNode=candidateNode,iNode=net[candidateNode].nn[candidateNodeIndex].nodeID;//determine link-end states
       	if(net[iNode].state==0)
       	{
		sNode=iNode;iNode=candidateNode;
	};

	net[sNode].state=1;					//SI->II
	varI++;							//update I-node count
	varSI+=net[sNode].nnState[0]-net[sNode].nnState[1];	//update SI-link count

        for(auto &i: net[sNode].nn)				//update neighborhoods of former S-node's neighbors
	{
        	net[i.nodeID].nnState[1]++;
        	net[i.nodeID].nnState[0]--;
	};
};

void rewiring(std::vector<netNode> &net,long int varI,long int &varSI,gsl_rng **rng,long int maxTries)//rewiring SI+S->SS+I
{	
	bool rewired=0;
	long int tries=0;
	long int n=net.size();
	while(rewired==0 && tries<maxTries && n-varI>1)//repeat if no previous success/not tried often enough/at least 1 S-node left
	{
		tries++;	//keep track of number of tries

//-------random pick of SI-link as in infection function 
		long int candidateLink=gsl_rng_uniform_int(*rng,2*varSI)+1;
		long int candidateNode=-1,siCounter=0;
		while(siCounter<candidateLink)
		{
			candidateNode++;
			if(net[candidateNode].state==0)
				siCounter+=net[candidateNode].nnState[1];
			else	
				siCounter+=net[candidateNode].nn.size()-net[candidateNode].nnState[1];
		};
		bool nodeState=net[candidateNode].state;
		long int candidateNodeIndex=net[candidateNode].nn.size()-1;
		while(siCounter>candidateLink || nodeState==net[net[candidateNode].nn[candidateNodeIndex].nodeID].state)
		{
			if(nodeState!=net[net[candidateNode].nn[candidateNodeIndex].nodeID].state)
				siCounter--;
			candidateNodeIndex--;
		};

		long int sNode=candidateNode,iNode=net[candidateNode].nn[candidateNodeIndex].nodeID;
       		if(net[iNode].state==0)
       		{
			sNode=iNode;iNode=candidateNode;
		};
//-----------------------------------------------------------------------

		long int sNode2=-1; 
		do					//randomly pick random S-node to rewire SI-link to
		{
			sNode2=-1;
			long int sCounter=0;
			long int pick=gsl_rng_uniform_int(*rng,n-varI)+1;
			while(sCounter<pick)
			{
				sNode2++;
				if(net[sNode2].state==0)
					sCounter++;
			};
		}
		while(sNode2==sNode);			//until picked S-node differs from S-end of picked SI-link
		bool linkExists=0;
		for(auto &i:net[sNode].nn)		//check if picked S-node is connected to S-end of picked SI-link
        	{
        		if(i.nodeID==sNode2) 
        	 	{
        	       		linkExists=1;
        	       		break;
        	        };
		};
        	if(linkExists==0)			//if all clear, do the actual rewiring
        	{

			removeLink(net,candidateNode,candidateNodeIndex);
        	        addLink(net,sNode,sNode2);
			varSI--;
        	        rewired=1;
        	        break;
      		};
	};
}
