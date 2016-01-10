//main.cpp is the control station - if functions.cpp is the API
//library, this is the place where you interact with it and
//call everything
#include <stdio.h> //simple output
#include <stdlib.h>
#include <mpi.h> //incorporates mpi capabilities
#include <fstream> //file writing
#include <iomanip> //setprecision()
#include <iostream>
#include <string>
#include <vector>
#include <omp.h> //extra parallelization
#include "bin_write.cpp"//read/write of binary files
#include "bin_read.cpp"
#include "functions.cpp"//all statistical functions, search class
//#include "structures.h"
//#include <gmp.h>
#define MVALS 8

//Import the list of log factorials
double *logFacts;
//Get the number of factorials, using the length
//of the binary list
int maxFact;

/*
double log_odds_ratio(double*,int,int*,int,double,double,bool); 
void normalize_counts(double*,int);
struct SearchResults;
void upload_data(double*,double*,int,int*,int*,int);
unsigned char *get_bins(double*,int,double*,int*,int*,int,double,double);
void free_data(double*,int*);
double bins_to_odds(unsigned char*,int,int*,int);
double t_odds(double*,int,double,double,int*,int);
double t_odds_two(double*,int,double,double,double,double,int,const char*);
*/

int main(int argc, char * argv[])
{
#pragma omp parallel for
    for (int i = 0; i < 100; i ++)
    {
        printf("Integer %d\n",i);
    }
	//read in some counts to practice on
	double *counts;
	//Get the total number of counts
	int length;
	//current search settings
	double nu, nudot;
	int m_max = MVALS;
	//settings/results packed into arrays:
	double curr_settings[2];
	double curr_results[3];
	//declare the mvalues to search
	int n_mvals = MVALS;
	int mvals[MVALS] = {2,4,8,16,32,64,128,255};
	//int mvals[8] = {2,4,8,16,32,64,128,255};
    printf("Run bushbaby, run!\n");
	//get filenames
    string filename = "data/data.bin";

    printf("Reading in logFacts\n");
	

	//load in values from data
	logFacts = bin_read((char*)"data/log_facs.bin");
	maxFact = bin_size((char*)"data/log_facs.bin");

	//counts = bin_read((char*)"data/B1821_counts.bin");
	//length = bin_size((char*)"data/B1821_counts.bin");
	//set up search
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
        printf("You entered %d args instead of 5 \n",argc-1);
		nu_min = 50; 
		nu_max = 500;
		nudot_min = 2.5e-19;//-1736.5e-16 
		nudot_max = 2.5e-11;//1736.5e-16;
	}
    counts = bin_read((char*)filename.c_str());
    length = bin_size((char*)filename.c_str());
    //normalize the counts
    normalize_counts(counts, length);
    string results_filename = filename.substr(0,filename.size()-3)+"results";
    printf("Process beginning search on %s\n"
            "Total %d counts spanning %f seconds\n"
            "Outputting best 5 odds to %s \n",
            filename.c_str(),
            length,counts[length-1],results_filename.c_str());


    int i = t_odds_two(counts, length, 
            nu_min,nu_max,
            nudot_min,nudot_max,
            verbosity,results_filename.c_str());
    if (!i)
    {
        printf("[+] Process completed search on %s"
                " with exit code %d\n", filename.c_str(),
                i);
    }
    else
    {
        printf("[-] Process exited search on %s"
                " with code %d\n", filename.c_str(),
                i);
    }
	printf("Process completed. Exiting.\n");
	free(logFacts);
	free(counts);
	return 0;
}
