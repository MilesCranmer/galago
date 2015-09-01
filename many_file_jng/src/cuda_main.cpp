//main.cpp is the control station - if functions.cpp is the API
//library, this is the place where you interact with it and
//call everything
#include <stdio.h> //simple output
#include <stdlib.h>
#include <gmp.h>
#include <mpi.h> //incorporates mpi capabilities
#include <fstream> //file writing
#include <iomanip> //setprecision()
#include <iostream>
#include <string>
#include <vector>
//#include <omp.h> //extra parallelization
#include "bin_write.cpp"//read/write of binary files
#include "bin_read.cpp"
//#include "functions.cu"//all statistical functions, search class
//#include "structures.h"
#define MVALS 8

//Import the list of log factorials
double *logFacts;
//Get the number of factorials, using the length
//of the binary list
int maxFact;

double log_odds_ratio(double*,int,int*,int,double,double,bool); 
void normalize_counts(double*,int);
struct SearchResults;
void upload_data(double*,double*,int,int*,int*,int);
unsigned char *get_bins(double*,int,double*,int*,int*,int,double,double);
void free_data(double*,int*);
double bins_to_odds(unsigned char*,int,int*,int);
double t_odds(double*,int,double,double,int*,int);
double t_odds_two(double*,int,double,double,double,double,int,const char*);

int main(int argc, char * argv[])
{
	//read in some counts to practice on
	double *counts;
	//Get the total number of counts
	int length;
	//current search settings
	double nu, nudot;
	int m_max = MVALS;
	//MPI variables (rank == proc number, size is num proc)
	int rank, size;
	//settings/results packed into arrays:
	double curr_settings[2];
	double curr_results[3];
	//declare the mvalues to search
	int n_mvals = MVALS;
	int mvals[MVALS] = {2,4,8,16,32,64,128,255};
	//int mvals[8] = {2,4,8,16,32,64,128,255};
	//holds all results (used by root node)
	//one communication channel - everything
	MPI_Comm comm = MPI_COMM_WORLD;
	//start up MPI
	MPI_Init(&argc, &argv);
	//get the process rank
	MPI_Comm_rank(comm,&rank);
	if (!rank)
		printf("Run baby run!\n");
	//get the number of processes
	MPI_Comm_size(comm,&size);
	//receiver status
	MPI_Status rstatus;
	MPI_Request srequest;
	//get filenames
	ifstream filenames_file("data/filenames.txt");
	vector<string> filenames;

	//char * srank;
	//sprintf(srank,"%d",rank);
	//string null_s(srank);
	//null_s += "_results.csv";
	//ofstream results_file(null_s.c_str());
	if (!rank)
		printf("Reading in files\n");

	string null_s = "";
	while (filenames_file.good())
	{
		for (int i = 0; i < rank; i++)
		{
			if (filenames_file.good())
				getline(filenames_file,null_s);
		}
		if (filenames_file.good())
			getline(filenames_file,null_s);
		if (null_s.length() == 0)
			break;
		filenames.push_back("data/"+null_s);
		for (int i = rank; i < size; i++)
		{
			if (filenames_file.good())
				getline(filenames_file,null_s);
		}
	}

	//as the root process, read in the data
	//if (rank == 0)
	if(!rank)
		printf("Reading in logFacts\n");
	

	//load in values from data
	logFacts = bin_read((char*)"data/log_facs.bin");
	maxFact = bin_size((char*)"data/log_facs.bin");

	//share sizing details - all of the following 
	//must be blocking to insure nothing funny happens
	/*
	if (size > 1)	
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
	*/
	//load in the data - these are blocking broadcasts, 
	//so act as barriers
	//if (size > 1)
		//MPI_Bcast(&logFacts[0], maxFact, MPI_DOUBLE, 0, comm);




	//split up MPI processes into different files
	//they pick the n*rank file, then the (n+1)*rank file,and so on.
	//counts = bin_read((char*)"data/B1821_counts.bin");
	//length = bin_size((char*)"data/B1821_counts.bin");
	//set up search
	if (!rank)
		printf("Applying search settings\n");
	int verbosity = 1;
	double nu_min,nu_max,nudot_min,nudot_max;
	if (argc == 6)
	{
		nu_min = atof(argv[1]);
		nu_max = atof(argv[2]);
		nudot_min = atof(argv[3]);//-1736.5e-16 
		nudot_max = atof(argv[4]);//1736.5e-16;
		verbosity = atoi(argv[5]);
	}
	else
	{
		if (!rank)
			printf("You entered %d args instead of 5 \n",argc-1);
		nu_min = 50; 
		nu_max = 500;
		nudot_min = 2.5e-19;//-1736.5e-16 
		nudot_max = 2.5e-11;//1736.5e-16;
	}
	if (!rank)
		printf("Parallel phase:\n");
	for (int file_i = 0; file_i < filenames.size(); file_i ++)
	{
		counts = bin_read((char*)filenames[file_i].c_str());
		length = bin_size((char*)filenames[file_i].c_str());
		//normalize the counts
		normalize_counts(counts, length);
		string results_filename = filenames[file_i].substr(0,filenames[file_i].size()-3)+"results";
		printf("Process %d beginning search on %s\n"
			   "Total %d counts spanning %f seconds\n"
			   "Outputting best 5 odds to %s \n",
			   rank,filenames[file_i].c_str(),
			   length,counts[length-1],results_filename.c_str());

		
		int i = t_odds_two(counts, length, 
						   nu_min,nu_max,
						   nudot_min,nudot_max,
						   verbosity,results_filename.c_str());
		if (!i)
		{
			printf("[+] Process %d completed search on %s"
				   " with exit code %d\n", rank,filenames[file_i].c_str(),
				   i);
		}
		else
		{
			printf("[-] Process %d exited search on %s"
				   " with code %d\n", rank,filenames[file_i].c_str(),
				   i);
		}

		//display some initial stats
		/*
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
		*/
		//go through all settings
		/*
		int i = 0;
		for (nu =  settings.nu_min;
			 nu <= settings.nu_max;
			 nu += settings.d_nu)
		{
			for (nudot =  settings.nudot_min;
				 nudot <= settings.nudot_max;
				 nudot += settings.d_nudot)
			{
				//double odds = log_odds_ratio(counts, length, m_max, 
										 	 //nu, nudot, 0);
				double odds = t_odds_two(counts, length,
										 nu, nudot);
				//some last modifications to the odds
				//finalodds = 1./nmvals/ log(nu_max/nu_min)/log(nudot_max/nudot_min) * dnu/nu * nudotfrac * moddsratio; 
				odds /= (settings.m_max - 1);
				odds *= settings.d_nu/nu;
				if (odds > 1e-8)
				{
					printf("Receiving Good Search, nu: %.9e, nudot: -%.9e, Odds: %.3e\n", 
					nu, nudot, odds);
					results_file << scientific << setprecision(10) << nu << ",";
					results_file << scientific << -nudot << ",";
					results_file << scientific << odds << ",";
					results_file << filenames[file_i] << ",";
					results_file << "\n";
					results.nu.push_back(curr_results[0]);
					results.nudot.push_back(curr_results[1]);
					results.odds.push_back(curr_results[2]);
				}
				i++;
			}
		}
		//end MPI
		results_file << "There were " << i << " total searches in " << filenames[file_i] << "\n";
		printf("Total searches: %d\n", i);
		*/
	}
	printf("Process %d completed. Exiting.\n", rank);
	MPI_Finalize();
	free(logFacts);
	free(counts);
	return 0;
}
