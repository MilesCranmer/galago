#include <math.h> //for math
#include <vector> //to hold search results
#include <algorithm> //compute max of vector
#include <numeric> //compute sum of vector (accumulate) 
//speed not important for final statistics, so optimising this is silly
#define PI 3.14159265359
//this file contains the main statistics functions
//This file relies on main's importing of the logFacts data file

//extern tells the compiler that this variable is global and
//is already initialized elsewhere
extern double *logFacts;
extern int maxFact;

using namespace std;

//Forward declaration for use in class
double log_odds_ratio(double *counts, int length, int m_max, 
					  double period, double phase, double nu, bool verbose);

//Holds the result of searches, and the displaying functions
//faster to have vectors of the variables within, 
//rather than a vector of the class
struct SearchResults
{
	//the settings of the search
	vector<double> period, phase;
	//The odds that this model fits
	vector<double> odds;
	//Get average and maximum odds
	double avg_odds();
	//return max odds
	double max_odds()
	{return *max_element(odds.begin(),odds.end());};
	//get the position of the max
	int max_odds_i()
	{return max_element(odds.begin(), odds.end()) - odds.begin();};
	//print settings for a specific search
	void print_settings(int i){printf("%lf\n",period[i]);};
	//print all settings
	void print_stats();
};

//This object organizes the entire search
struct PeakSearch
{
	//the bounds on the period search
	double period_min,period_max; 
	//interval between searches w.r.t period
	double d_period; 
	//the max number of bins
	int m_max; 
	//the bounds on the phase search
	//**In terms of radians
	double phase_min, phase_max;
	//interval between searches w.r.t phase
	double d_phase; 
	//holds the value of nu, which I still don't 
	//really understand (hence I set to 1)
	double nu;
	//set default parameters
	void default_params();
	//Print all of the above parameters
	void print_params();
	SearchResults search(double *counts, int length, bool verbose);
};

//gets the average odds
double SearchResults::avg_odds()
{
	//compute the sum of the odds ratios
	double sum = accumulate(odds.begin(), odds.end(),0);
	//get the average, and return
	return sum/odds.size();
}

//print some stats about the search
void SearchResults::print_stats()
{
	printf("\nSTATS:\n");
	printf("The average odds for the search ");
	printf("are %lf\n",avg_odds());
	//Get position of max
	int max_i = max_odds_i();
	printf("The best odds occur for ");
	printf("a period of %lf seconds ",period[max_i]);
	printf("and a phase of %lf radians, which ",phase[max_i]);
	printf("give odds of %lf\n",odds[max_i]);
}

//set defaults
void PeakSearch::default_params()
{
	//the bounds on the period search
	period_min = 5;
	period_max = 10;
	//interval between searches w.r.t period
	d_period = (period_max-period_min)/30;
	//the max number of bins
	m_max = 50;
	//the bounds on the phase search
	//**In terms of radians, hence
	//0->2 pi is a full search
	phase_min = 0;
	//interval between searches w.r.t phase
	d_phase = PI/3;
	//don't repeat first phase
	phase_max = 2*PI-d_phase;
	//set to value that won't effect outcome
	nu = 1;
}

//Print out all settings
void PeakSearch::print_params()
{
	printf("The period is tested from %lf to %lf seconds\n",
		   period_min,period_max);
	printf("The interval of this search is %lf seconds\n",d_period);
	printf("The phase is tested from %lf to %lf radians\n",
		   phase_min,phase_max);
	printf("The interval of this search is %lf radians\n",d_phase);
	printf("The maximum number of bins in the stepwise model is %d\n"
		,m_max);
}
//double log_odds_ratio(double *counts, int length, int m_max, 
//					  double period, double phase, double nu)
//searches through all settings
//length is the number of counts, and verbose
//tells the program how much it should talk
SearchResults PeakSearch::search(double *counts, int length,
								 bool verbose)
{
	//holds all search settings
	//an array of three element arrays.
	SearchResults searches;
	if(verbose)
		printf("Starting searches\n");
	//iterate through possible periods and phases
	for (double period = period_min; period <= period_max; 
		 period += d_period)
	{
		for (double phase = phase_min; phase <= phase_max;
			 phase += d_phase)
		{
			//each setting will be in the same index,
			//so can be accessed later
			searches.period.push_back(period);
			searches.phase.push_back(phase);
			double odds = log_odds_ratio(counts, length, m_max, 
										 period, phase, nu, verbose);
			searches.odds.push_back(odds);
			if(verbose)
			{	
				printf("period=%lf,phase=%lf,odds=%lf\n",
						period,phase,odds);
			}
		}
	}
	//return all computed searches
	return searches;
}

//This function returns the choose function
double log_choose(int first, int second)
{
	//sanity check
	if (first > maxFact || second > maxFact || second > first)
	{
		return 0;
	}
	//the log of the choose function
	return logFacts[first] - logFacts[second] - logFacts[first-second];
}

//this function normalizes a list of counts to start at 0 and 
//be in terms of seconds
void normalize_counts(double *counts, int length)
{
	double t_min = counts[0];
	for (int i = 0; i < length; i ++)
	{
		counts[i] -= t_min;
	}
}

//This function gets the number of events which occur within
//the specified range. Length is the total number of events
int num_events(double *counts, int length, double start, double end)
{
	//the number of events in this region of the counts
	int num = 0;
	//go through all counts
	//note that this assumes they are in order
	for (int i = 0; i < length; i ++)
	{
		//check if the count is in the region
		if (counts[i] >= start && counts[i] < end)
		{
			num ++;
		}
		//if the count is later, then we choose to not be redundant
		else if (counts[i] >= end)
		{
			break;
		}
	}
	//return the number of events
	return num;
}

//Equation from gregory and loredo paper to calcluate odds ratio
//of m-binned stepwise model w.r.t. constant model
double log_m_odds_ratio(double *counts, int length, int m, 
					  double period, double phase_eff, double nu,
					  double t_max)
{
	//They will be used to normalize all the count times
	double binWidth = period/((double)m);
	//All of the nj log factorials together
	double allnj = 0;
	//the final log ratio
	double om1 = 0;

	//go through all bins
	for (int j = 0; j < m; j++)
	{
		//number of counts in this bin
		int nj = 0;
		//start of the first portion of this bin
		double start = binWidth*((double)j)-period+phase_eff;
		//Go through all portions of this bin (folded period)
		for (double t_curr = start; t_curr <= t_max+period+phase_eff; 
			t_curr += period)
		{
			//get the number of events within this portion
			nj += num_events(counts, length, t_curr, 
							 t_curr + binWidth);
		}
		allnj += logFacts[nj];

	}
	//calculate the final log ratio
	om1 += allnj;
	om1 += ((double)length)*log(m);
	om1 -= log_choose(length+m-1,length);
	om1 -= logFacts[length] + log(nu);
	return om1;
}

//Equation from gregory and loredo paper to calcluate total odds
//ratio
double log_odds_ratio(double *counts, int length, int m_max, 
					  double period, double phase, double nu, bool verbose)
{
	//normalize the counts with item 0 at t=0.0s
	normalize_counts(counts, length);
	//the following assumes the counts are ordered
	double t_max = counts[length-1];
	//The total odds ratio
	double odds = 0;
	//the actual phase in seconds
	double phase_eff = phase/2/PI*period;
	//go through all possible m values
	for (int m = 2; m <= m_max; m++)
	{
		if (verbose)
			printf("Testing %d-binned model\n",m);
		//Add the next om1 value to the total odds ratio.
		//We also have to remove the log
		odds += exp(log_m_odds_ratio(counts,length,m,period,
									 phase_eff,nu,t_max));
	}
	return odds;
}

//Gets the average time between counts
double avg_interval(double *counts, int length)
{
	double total_time;
	total_time = counts[length-1] - counts[0];
	return total_time/length;
}

//Gets the minimum time between counts
double min_interval(double *counts, int length)
{
	double smallest;
	//start with first interval
	smallest = counts[1] - counts[0];
	//go through the rest
	for (int i = 2; i < length; i ++)
	{
		double tmp = counts[i] - counts[i-1];
		//if interval smaller, assume it is now
		//the smallest
		if (tmp < smallest)
		{
			smallest = tmp;
		}
	}
	return smallest;
}