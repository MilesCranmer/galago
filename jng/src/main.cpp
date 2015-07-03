//main.cpp is the control station - if functions.cpp is the API
//library, this is the place where you interact with it and
//call everything
#include <stdio.h> //simple output
#include <stdlib.h>
#include <mpi.h> //incorporates mpi capabilities
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
	double period, phase;
	double nu = 1;
	int m_max = 50;
	//MPI variables (rank == proc number, size is num proc)
	int rank, size;
	//settings/results packed into arrays:
	double curr_settings[2];
	double curr_results[3];
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
		logFacts = bin_read((char*)"data/log_facs.bin");
		maxFact = bin_size((char*)"data/log_facs.bin");
		counts = bin_read((char*)"data/practice_counts.bin");
		length = bin_size((char*)"data/practice_counts.bin");
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
	//printf("Proc %d: Log(50!) = %f, num of facts = %d, count 100 = %f\n", 
	//		rank, logFacts[50], maxFact, counts[99]);


	if (rank == 0)
	{
		//holds all results
		SearchResults results;
		//set up search
		PeakSearch settings;
		//call default parameters
		settings.default_params();
		settings.d_period = 1.00;

		//display some initial stats
		printf("The settings for this search are:\n\n");
		printf("Min Period: %f s, Max Period: %f s\n",
			   settings.period_min, settings.period_max);
		printf("Period Interval: %f s\n",
			   settings.d_period);
		printf("Min Phase: %f rad, Max Phase: %f rad\n",
			   settings.phase_min, settings.phase_max);
		printf("Phase Interval: %f rad\n", 
			   settings.d_phase);
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
		for (period =  settings.period_min;
			 period <= settings.period_max;
			 period += settings.d_period)
		{
			//new settings to send
			curr_settings[0] = period;
			for (phase =  settings.phase_min;
				 phase <= settings.phase_max;
				 phase += settings.d_phase)
			{
				//new settings to send
				curr_settings[1] = phase;
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
					MPI_Recv(curr_results, 3, MPI_DOUBLE, i, 
							 chan_results, comm, &rstatus);
					//one more search received
					recv ++;
					printf("Receiving Search, Period: %f, Phase: %f, Odds: %f\n", 
							curr_results[0], curr_results[1], curr_results[2]);
					results.period.push_back(curr_results[0]);
					results.phase.push_back(curr_results[1]);
					results.odds.push_back(curr_results[2]);
				}
				//printf("Sending Search to proc %d, Period: %f, Phase: %f\n", 
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
			MPI_Recv(curr_results, 3, MPI_DOUBLE, i, chan_results, comm, &rstatus);
			//print results
			printf("Receiving Search, Period: %f, Phase: %f, Odds: %f\n", 
					curr_results[0], curr_results[1], curr_results[2]);
			results.period.push_back(curr_results[0]);
			results.phase.push_back(curr_results[1]);
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
		i = results.max_odds_i();
		printf("\nThe best odds are: %f, which occur for a period of %f and"
			   " phase of %f\n\n",
				results.odds[i], results.period[i], results.phase[i]);
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
										 	nu, 0);
				curr_results[0] = curr_settings[0];
				curr_results[1] = curr_settings[1];
				//send back results of search
				MPI_Isend(curr_results, 3, MPI_DOUBLE,
						  0, chan_results, comm, &srequest);
				data_flag = 0;
			}
			//check for shutdown
			MPI_Iprobe(0, chan_terminate, comm, &termination_flag, &rstatus);
		}
	}
	/*
	PeakSearch settings;
	settings.default_params();
	settings.period_min = 0.001;
	settings.period_max = 0.010;
	settings.d_period = 0.001;
	settings.m_max = 20;
	settings.print_params();
	SearchResults searches;
	searches = settings.search(counts, length, true);
	searches.print_stats();*/

	printf("proc %d terminating\n", rank);
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