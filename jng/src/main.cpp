//main.cpp is the control station - if functions.cpp is the API
//library, this is the place where you interact with it and
//call everything
#include <stdio.h> //simple output
#include <stdlib.h>
#include <gmp.h>
#include <mpi.h> //incorporates mpi capabilities
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
	//read in some counts to practice on
	double *counts;
	//Get the total number of counts
	int length;
	//current search settings
	double nu, nudot;
	int m_max = 100;
	//MPI variables (rank == proc number, size is num proc)
	int rank, size;
	//settings/results packed into arrays:
	double curr_settings[2];
	double curr_results[3];
	//holds all results (used by root node)
	SearchResults results;
	//whether or not start of 
	bool initialized = 0;
	//flags used in communication
	//data checks if new data available, while
	//termination flag tells processes to terminate
	int data_flag = 0, termination_flag = 0;
	//two channels are for data to be sent, the other
	//is to check if terminated
	int chan_results = 2, chan_settings = 1;
	int chan_terminate = 0;
	//one communication channel - everything
	MPI_Comm comm = MPI_COMM_WORLD;
	//start up MPI
	MPI_Init(&argc, &argv);
	//get the process rank
	MPI_Comm_rank(comm,&rank);
	//get the number of processes
	MPI_Comm_size(comm,&size);
	//receiver status
	MPI_Status rstatus;
	MPI_Request srequest;

	//as the root process, read in the data
	if (rank == 0)
	{
		//load in values from data
		logFacts = bin_read((char*)"data/log_facs_2.bin");
		maxFact = bin_size((char*)"data/log_facs_2.bin");
		counts = bin_read((char*)"data/B1821_counts.bin");
		length = bin_size((char*)"data/B1821_counts.bin");
		//normalize the counts
		normalize_counts(counts, length);
		printf("Total of %d counts\n", length);
		//print last count
		printf("After normalization, the last count is at %lf seconds\n\n", counts[length-1]);
	}

	//share sizing details - all of the following 
	//must be blocking to insure nothing funny happens
	MPI_Bcast(&maxFact, 1, MPI_INT, 0, comm);
	MPI_Bcast(&length, 1, MPI_INT, 0, comm);
	if(rank != 0)
	{
		//allocate memory to hold the entered data
		/////////////////////////////////////////////////
		//If a guillimin node cannot hold all of the toas,
		//this part will instead be transferring the data
		//to a binary file, which can be read over and over
		/////////////////////////////////////////////////
		logFacts = new double[maxFact];
		counts = new double[length];
	}
	//load in the data - these are blocking broadcasts, 
	//so act as barriers
	MPI_Bcast(&logFacts[0], maxFact, MPI_DOUBLE, 0, comm);
	MPI_Bcast(&counts[0], length, MPI_DOUBLE, 0, comm);
	//make sure all processes have the correct data
	//printf("Proc %d: Log(50!) = %e, num of facts = %d, count 100 = %e\n", 
	//		rank, logFacts[50], maxFact, counts[99]);


	if (rank == 0)
	{
		//set up search
		PeakSearch settings;
		//call default parameters
		settings.default_params();
		settings.nu_min = 325;
		settings.nu_max = 330;
		settings.d_nu = 0.00001;
		settings.nudot_min = 1.7365e-15; 
		settings.nudot_max = 1.7365e-15;
		settings.d_nudot = 1e-16;
		settings.m_max = 50;

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

		//iterator to go through all processes in an
		//even fashion (equal distribution of work)
		int i = 1;
		int sent = 0;
		int recv = 0;
		//go through all settings
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
				//new settings to send
				curr_settings[1] = nudot;
				//once everything is up and running, go through
				//work loop
				if (initialized)
				{
					//continually check all processes
					//until one sends data.
					while(!data_flag)
					{
						//check if program is finished by seeing
						//if results have been sent yet
						MPI_Iprobe(i, chan_results, comm, &data_flag, 
								   &rstatus);
						//if still no data
						if (!data_flag)
						{
							//go to next proc
							i++;
							//modulate the iterator
							if(i == size){i = 1;}
						}
					}
					//data was sent to root:
					//Collect search results (blocking as we already know 
					//a message was sent from probe)
					MPI_Recv(curr_results, 3, MPI_DOUBLE, i, 
							 chan_results, comm, &rstatus);

					//one more search received
					recv ++;
					//some last modifications to the odds
					//finalodds = 1./nmvals/ log(nu_max/nu_min)/log(nudot_max/nudot_min) * dnu/nu * nudotfrac * moddsratio; 
					curr_results[2] /= (settings.m_max - 1);
					curr_results[2] *= settings.d_nu/nu;
					//printf("Curr %e\n", curr_results[2]);
					//curr_results[2] /= log(settings.nu_max/settings.nu_min);
					//printf("Curr %e\n", curr_results[2]);
					//curr_results[2] *= settings.d_nu/nu;
					//printf("Curr %e\n", curr_results[2]);
					//curr_results[2] *= fabs(settings.d_nudot/nudot);
					printf("Receiving Search, nu: %e, nudot: -%e, Odds: %e\n", 
							curr_results[0], curr_results[1], curr_results[2]);
					results.nu.push_back(curr_results[0]);
					results.nudot.push_back(curr_results[1]);
					results.odds.push_back(curr_results[2]);
				}
				//printf("Sending Search to proc %d, nu: %e, nudot: %e\n", 
				//	   i, curr_settings[0], curr_settings[1]);
				//send next settings (non-blocking so can keep working)
				MPI_Isend(curr_settings, 2, MPI_DOUBLE, i, 
						  chan_settings, comm, &srequest);
				//one more search started
				sent ++;

				//awaiting next data
				data_flag = 0;
				//go to next proc
				i++;
				//check if all processes started, and if true, 
				//then the program has been initialized
				if (!initialized && i == size){initialized = 1;}
				//modulate the iterator
				if(i == size){i = 1;}
			}
		}
		//collect rest of results
		while (recv < sent)
		{
			//continually check all processes
			//until one sends data.
			while(!data_flag)
			{
				//check if program is finished by seeing
				//if results have been sent yet
				MPI_Iprobe(i, chan_results, comm, &data_flag, &rstatus);
				//if still no data
				if (!data_flag)
				{
					//go to next proc
					i++;
					//modulate the iterator
					if(i == size){i = 1;}
				}
			}
			//data was sent to root: 
			//Collect search results (blocking as we already know 
			//a message was sent from probe)
			MPI_Recv(curr_results, 3, MPI_DOUBLE, i, chan_results, comm, 
					 &rstatus);
			//print results
			printf("Receiving Search, nu: %e, nudot: -%e, Odds: %e\n", 
					curr_results[0], curr_results[1], curr_results[2]);
			results.nu.push_back(curr_results[0]);
			results.nudot.push_back(curr_results[1]);
			results.odds.push_back(curr_results[2]);
			//one more search received
			recv ++;
			//set back data_flag
			data_flag = 0;
		}
		//end everything
		for (i = 1; i < size; i++)
		{
			//variable sent does not matter
			//this will trigger
			//termination of all other processes
			//than root
			MPI_Isend(&size, 1, MPI_INT,
					  i, chan_terminate, comm,
					  &srequest);
		}
	}
	else //the bulk search processes
	{
		//go until told to terminate
		while(!termination_flag)
		{
			//check for new settings
			MPI_Iprobe(0, chan_settings, comm, &data_flag, &rstatus);
			if (data_flag)
			{
				//first, receive settings and unpack
				MPI_Recv(curr_settings, 2, MPI_DOUBLE,
						 0, chan_settings, comm, &rstatus);
				//new search settings. Perform search.
				///////////////////////////////////////
				//Eventually OpenMP should do something
				//to parallelize the log odds ratio
				//////////////////////////////////////s
				curr_results[2] = log_odds_ratio(counts, length, m_max, 
										 	curr_settings[0], 
										 	curr_settings[1], 
										 	0);
				curr_results[0] = curr_settings[0];
				curr_results[1] = curr_settings[1];
				//send back results of search
				MPI_Isend(curr_results, 3, MPI_DOUBLE,
						  0, chan_results, comm, &srequest);
				data_flag = 0;
			}
			//check for shutdown
			MPI_Iprobe(0, chan_terminate, comm, &termination_flag, 
					   &rstatus);
		}
	}
	if (rank != 0)
	{
		printf("proc %d terminating\n", rank);
	}
	//so results message is last
	MPI_Barrier(comm);
	if (rank == 0)
	{
		int i = results.max_odds_i();
		printf("\nThe best odds are: %e, which occur for nu of %e Hz and"
			   " nudot of -%e Hz/s\n\n",
				results.odds[i], results.nu[i], results.nudot[i]);
	}
	//end MPI
	MPI_Finalize();
	if (rank != 0)
	{
		//free the worker process memory of
		//the logarithmic factorials and list of counts
		free(logFacts);
		free(counts);
	}
	return 0;
}
