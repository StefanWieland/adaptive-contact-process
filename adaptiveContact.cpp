#include "adaptiveContactFunctions.h"

int main(int ac,char **av)
{
	time_t t1,t2;(void)time(&t1);			//initializing time
//------------------reading in parameters---------------------
	const long int N=atol(av[1]);		//number of network nodes
	const long int K=atol(av[2]);		//number of network links
	const long int I0=atol(av[3]);		//number of initially infected nodes
	const double P=atof(av[4]);		//transmission rate SI->II
	const double R=atof(av[5]);		//recovery rate I->S
	const double W=atof(av[6]);		//rewiring rate SI+S->SS+I
	const long int TMAX=atol(av[7]);	//simulated time
	const long int RUNS=atol(av[8]);	//number of runs to average over
	const long int TRIESMAX=atol(av[9]);	//maximum number of tries to carry out rewiring event
	const long int KMAX=atol(av[10]);	//maximum considered degree
	std::string type=av[11];		//file name 
//-------------------------------------------------------------

	long int varIglobal=0;			//counter for averaging [I] (densities of I-nodes)
	long int varSIglobal=0;			//counter for averaging [SI] (densitites of SI-links)
	double Deg[KMAX][KMAX][2]={};		//joint degree distribution Deg[S-neighbors][I-neighbors][P_S,P_I] 

//------------------prepare status bar visualization----------------
	long int runCounter=0, starCounter=0;	
	std::cout<<"\n0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%	progress\n";	
	std::cout<<"|    |    |    |    |    |    |    |    |    |    |\n";
	std::cout<<"*";
//-----------------------------------------------------------------

	double t0=time(0);			//time stamp;different rng seeds for different program executions 

	#pragma omp parallel for reduction(+:varIglobal,varSIglobal)//
	for(long int run=0;run<RUNS;++run)			//averaging over runs
	{
		long int varI=0,varSI=0;			//propensity counters

		gsl_rng *rng=gsl_rng_alloc(gsl_rng_default);	//declaring thread-specific rng
		gsl_rng_set(rng,t0+run);			//thread-safe seeding of rng
	
		std::vector<netNode> net(N);			//declaring network vector
		stateInitialization(net,I0,&rng);		//initialize with I0 I-nods	
		ErdosRenyi(net,K,&rng);				//choose Erdos-Renyi graph
		countMoments(net,varI,varSI);			//count necessary densities for Gillespie algorithm
//		printNetwork(net);				//print network

		double t=0;					//initialize system time		

		//simulate stochastic trajectory with Gillespie algorithm.stop if time's up,no I-nodes or no SI-links left.
		while(t<TMAX && varI>0 && varSI>0)
		{
			double process=gsl_rng_uniform(rng)*(R*varI+(W+P)*varSI); //weight

			double dtrand=0;					  //time step
			while(dtrand==0)
				dtrand=gsl_rng_uniform(rng);
			double dt=-1/(R*varI+(W+P)*varSI)*log(dtrand);
       		 	t+=dt;

			if(process<R*varI)					  //choosing reaction channels
	   		    	 recovery(net,varI,varSI,&rng);			
	    		else
				if(process<R*varI+P*varSI)
					infection(net,varI,varSI,&rng);
				else
					if(W!=0)
						rewiring(net,varI,varSI,&rng,TRIESMAX);
		};
		#pragma omp atomic	//updating density counters
		varIglobal+=varI;
		#pragma omp atomic
		varSIglobal+=varSI;

		for(auto &i : net) 	//record joint degree distributions
		{		
			long int x=i.nn.size()-i.nnState[1], y=i.nnState[1];//always >=0
			int status=i.state;
			if(x<KMAX && y<KMAX)
			{
				#pragma omp atomic
				Deg[x][y][status]++;
			};
		};
		(void)time(&t2);

//-----------status bar visualization---------------------------------------
		#pragma omp atomic
		runCounter++;
		#pragma omp critical
		{
			if(RUNS>=50)
			{
				if(runCounter>=double(RUNS*(starCounter+1))/50)
				{
					std::cout<<"*"<<std::flush;;
					starCounter++;
				};
			}
			else
			{
				int mustHaveStars=int(50*runCounter/RUNS);
				int addStars=mustHaveStars-starCounter;
				starCounter=mustHaveStars;
				for(int j=0;j<addStars;++j)
					std::cout<<"*"<<std::flush;;
			};
		};
//-------------------------------------------------------------------------------

		gsl_rng_free(rng);	//freeing allocated memory of thread-specific rng		
	};

//----writing joint degree distributions----------------------
	std::vector<double> norm(2);				//computing normalizations
	for(long int x=0;x<KMAX;++x)
		for(long int y=0;y<KMAX;++y)
			for(int j=0;j<2;++j)
				norm[j]+=Deg[x][y][j];

	(void)time(&t2);					//2nd flag of computation time 
	std::fstream data;
	data.open((type).c_str(), std::ios::out);		//general data
	data<<"# average [I]: "<<double(varIglobal)/N/RUNS<<", average [SI]: "<<double(varSIglobal)/N/RUNS<<", computation time: "<<t2-t1<<"s\n";
	for(long int x=0;x<KMAX;++x)					//browse S-neighbors
	{
		for(long int y=0;y<KMAX;++y)				//browse I-neighbors
			data<<x<<" "<<y<<" "<<Deg[x][y][0]/norm[0]<<" "<<Deg[x][y][1]/norm[1]<<"\n";
		data<<"\n";
	};
	data.close();
	std::cout<<"\n\n";
}


