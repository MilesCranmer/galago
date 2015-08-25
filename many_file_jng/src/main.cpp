//main.cpp is the control station - if functions.cpp is the API
//library, this is the place where you interact with it and
//call everything
#include <stdio.h> //simple output
#include <stdlib.h>
#include <gmp.h>
#include <mpi.h> //incorporates mpi capabilities
#include <fstream> //file writing
#include <iomanip> //setprecision()
//#include <omp.h> //extra parallelization
#include "bin_write.cpp"//read/write of binary files
#include "bin_read.cpp"
#include "functions.cpp"//all statistical functions, search class
//Import the list of log factorials
double *logFacts;
//Get the number of factorials, using the length
//of the binary list
int maxFact;

//int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
//             MPI_Comm comm)
//int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
//             MPI_Comm comm, MPI_Status *status)

int main(int argc, char * argv[])
{
	//results file
	ofstream results_file("del_me.csv");
	//read in some counts to practice on
	double *counts;
	//Get the total number of counts
	int length;
	//current search settings
	double nu, nudot;
	int m_max = 15;
	//MPI variables (rank == proc number, size is num proc)
	int rank, size;
	//settings/results packed into arrays:
	double curr_settings[2];
	double curr_results[3];
	//holds all results (used by root node)
	SearchResults results;
	//one communication channel - everything
	MPI_Comm comm = MPI_COMM_WORLD;
	//start up MPI
	MPI_Init(&argc, &argv);
	//get the process rank
	MPI_Comm_rank(comm,&rank);
	//get the number of processes
	MPI_Comm_size(comm,&size);

	//as the root process, read in the data
	if (rank == 0)
	{
		//load in values from data
		logFacts = bin_read((char*)"data/log_facs_2.bin");
		maxFact = bin_size((char*)"data/log_facs_2.bin");
	}
	//share sizing details - all of the following 
	//must be blocking to insure nothing funny happens
	MPI_Bcast(&maxFact, 1, MPI_INT, 0, comm);
	if(rank != 0)
	{
		//allocate memory to hold the entered data
		/////////////////////////////////////////////////
		//If a guillimin node cannot hold all of the toas,
		//this part will instead be transferring the data
		//to a binary file, which can be read over and over
		/////////////////////////////////////////////////
		logFacts = new double[maxFact];
	}
	//load in the data - these are blocking broadcasts, 
	//so act as barriers
	MPI_Bcast(&logFacts[0], maxFact, MPI_DOUBLE, 0, comm);
	//make sure all processes have the correct data
	//printf("Proc %d: Log(50!) = %e, num of facts = %d, count 100 = %e\n", 
	//		rank, logFacts[50], maxFact, counts[99]);

	//split up MPI processes into different files
	//they pick the n*rank file, then the (n+1)*rank file,and so on.
	if (rank == 0)
	{
		counts = bin_read((char*)"data/B1821_counts.bin");
		length = bin_size((char*)"data/B1821_counts.bin");
	}
	else
	{
		free(logFacts);
		MPI_Finalize();
		return 0;
	}
	//normalize the counts
	normalize_counts(counts, length);


	//set up search
	PeakSearch settings;
	//call default parameters
	settings.default_params();
	settings.nu_min = 327.2;
	settings.nu_max = 327.5;
	settings.d_nu = 1/counts[length-1];
	settings.nudot_min = 1736.5e-16;//-1736.5e-16 
	settings.nudot_max = 1736.5e-16;
	settings.d_nudot = 1e-8;
	settings.m_max = 15;

	//display some initial stats
	printf("The settings for this search are:\n\n");
	printf("Min nu: %e Hz, Max nu: %e Hz\n",
		   settings.nu_min, settings.nu_max);
	printf("nu Interval: %e Hz\n",
		   settings.d_nu);
	printf("Min nudot: %e Hz/s, Max nudot: %e Hz/s\n",
		   settings.nudot_min, settings.nudot_max);
	printf("nudot Interval: %e Hz/s\n", 
		   settings.d_nudot);
	printf("Max number of bins: %d\n", 
		   settings.m_max);

	printf("\n*******************************\n");
	printf("Starting search!\n");
	printf("*******************************\n\n");
	//go through all settings
	int i = 0;
	for (nu =  settings.nu_min;
		 nu <= settings.nu_max;
		 nu += settings.d_nu)
	{
		//new settings to send
		curr_settings[0] = nu;
		for (nudot =  settings.nudot_min;
			 nudot <= settings.nudot_max;
			 nudot += settings.d_nudot)
		{
			double odds = log_odds_ratio(counts, length, m_max, 
									 	 nu, nudot, 0);
			//some last modifications to the odds
			//finalodds = 1./nmvals/ log(nu_max/nu_min)/log(nudot_max/nudot_min) * dnu/nu * nudotfrac * moddsratio; 
			odds /= (settings.m_max - 1);
			odds *= settings.d_nu/nu;
			if (odds > 1e-2)
			{
				printf("Receiving Good Search, nu: %.9e, nudot: -%.9e, Odds: %.3e\n", 
				nu, nudot, odds);
				results_file << scientific << setprecision(10) << nu << ",";
				results_file << scientific << -nudot << ",";
				results_file << scientific << odds << ",";
				results_file << "\n";
				results.nu.push_back(curr_results[0]);
				results.nudot.push_back(curr_results[1]);
				results.odds.push_back(curr_results[2]);
			}
			i++;
		}
	}
	//end MPI
	results_file << "There were " << i << " total searches\n";
	printf("Total searches: %d\n", i);
	MPI_Finalize();
	free(logFacts);
	free(counts);
	return 0;
}