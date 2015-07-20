#include <math.h> //for math
//#include <gmpxx.h> //for precision calculation
#include <vector> //to hold search results
#include <algorithm> //compute max of vector
#include <numeric> //compute sum of vector (accumulate)
#ifdef _OPENMP
#include <omp.h>
#endif
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
					  double nu, double nudot, bool verbose);

//Holds the result of searches, and the displaying functions
//faster to have vectors of the variables within, 
//rather than a vector of the class
struct SearchResults
{
	//the settings of the search
	vector<double> nu, nudot;
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
	void print_settings(int i){printf("%lf\n",nu[i]);};
	//print all settings
	void print_stats();
};

//This object organizes the entire search
struct PeakSearch
{
	//the bounds on the nu search
	double nu_min,nu_max; 
	//interval between searches w.r.t nu
	double d_nu; 
	//the max number of bins
	int m_max; 
	//the bounds on the nudot search
	//**In terms of radians
	double nudot_min, nudot_max;
	//interval between searches w.r.t nudot
	double d_nudot; 
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
	printf("a nu of %lf seconds ",nu[max_i]);
	printf("and a nudot of %lf radians, which ",nudot[max_i]);
	printf("give odds of %lf\n",odds[max_i]);
}

//set defaults
void PeakSearch::default_params()
{
	//the bounds on the nu search
	//in terms of Hz
	nu_min = 5;
	nu_max = 10;
	//interval between searches w.r.t nu
	d_nu = (nu_max-nu_min)/30;
	//the max number of bins
	m_max = 30;
	//the bounds on the nudot search
	//In terms of Hz/s
	nudot_min = 1e-100;
	//don't repeat first nudot
	nudot_max = 2e-100;
	//interval between searches w.r.t nudot
	d_nudot = 1e-100;
}

//Print out all settings
void PeakSearch::print_params()
{
	printf("The nu is tested from %lf to %lf seconds\n",
		   nu_min,nu_max);
	printf("The interval of this search is %lf seconds\n",d_nu);
	printf("The nudot is tested from %lf to %lf radians\n",
		   nudot_min,nudot_max);
	printf("The interval of this search is %lf radians\n",d_nudot);
	printf("The maximum number of bins in the stepwise model is %d\n"
		,m_max);
}
//double log_odds_ratio(double *counts, int length, int m_max, 
//					  double nu, double nudot, double nu)
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
	//iterate through possible nus and nudots
	for (double nu = nu_min; nu <= nu_max; 
		 nu += d_nu)
	{
		for (double nudot = nudot_min; nudot <= nudot_max;
			 nudot += d_nudot)
		{
			//each setting will be in the same index,
			//so can be accessed later
			searches.nu.push_back(nu);
			searches.nudot.push_back(nudot);
			//finalodds = 1./nmvals/ log(nu_max/nu_min)/log(nudot_max/nudot_min) * dnu/nu * nudotfrac * moddsratio; 
			double odds = 1;
			odds *= 1/m_max/log(nu_max/nu_min)/log(nudot_max/nudot_min)*d_nu/nu*fabs(d_nudot/nudot);
			odds *= log_odds_ratio(counts, length, m_max, 
								   nu, nudot, verbose);
			printf("Odds: %f\n", odds);
			searches.odds.push_back(odds);
			if(verbose)
			{	
				printf("nu=%lf,nudot=%lf,odds=%lf\n",
						nu,nudot,odds);
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
					  double nu, double nudot,
					  double t_max)
{
	//odds to return
	double om1 = 0.0;
	//create all the bins, init to zero counts
	int n[m];

	/*
	//declare variables (use of which is declared later on)
	mpf_t  mpg_time, mpg_phi, mpg_nu, mpg_nudot, mpg_trunc_phi, mpg_d_phi; 
	//initialize them!
	mpf_init(mpg_time); 
	mpf_init(mpg_phi); 
	mpf_init(mpg_nu); 
	mpf_init(mpg_nudot); 
	mpf_init(mpg_trunc_phi); 
	mpf_init(mpg_d_phi); 

	//set search values
	mpf_set_d(mpg_nu, nu);
	mpf_set_d(mpg_nudot, nudot);
	*/

	//init to zero
	for (int j = 0; j < m; j++)
	{
		n[j] = 0;
	}

	//variables used in binnings
	//gets position in nu
	long double phi, d_phi;
	//double phi;
	//gets bin
	int k;
	//bin the photons
#ifdef _OPENMP
#pragma omp parallel for private(phi, d_phi, k) shared(n)
#endif
	for (int i = 0; i < length; i++)
	{
		
		/*
		//set current time
		mpf_set_d(mpg_time, counts[i]);
		//calculate changing of bin (dphi = 0.5*t^2*nudot)
		mpf_mul(mpg_d_phi, mpg_time, mpg_time);
		mpf_mul(mpg_d_phi, mpg_d_phi, mpg_nudot);
		mpf_div_2exp(mpg_d_phi, mpg_d_phi, 1);
		//calculate phase of bin, using changing of bin
		//(phi = (t*nudot+dphi)mod1)
		mpf_mul(mpg_phi, mpg_time, mpg_nu);
		mpf_add(mpg_phi, mpg_phi, mpg_d_phi);
		//do the mod1 part
		mpf_trunc(mpg_trunc_phi, mpg_phi);
		mpf_sub(mpg_phi, mpg_phi, mpg_trunc_phi);
		phi = mpf_get_d(mpg_phi);
		//get the corresponding bin
		k = (int)(phi*m);
		//one more count to the bin!
		n[k]++;
		*/
		
		//old normal precision calculation
		//in period
		
		d_phi = 0.5*counts[i]*nudot*counts[i];
		//get position in nu of photon
		phi = fmod(counts[i]*nu+d_phi,1);
		//get corresponding bin	
		k = (int)(phi*m);
		//one more count
		n[k]++;
		
	}

	/*
	//clean up everything
	mpf_clear(mpg_time); 
	mpf_clear(mpg_phi); 
	mpf_clear(mpg_nu); 
	mpf_clear(mpg_nudot); 
	mpf_clear(mpg_trunc_phi); 
	mpf_clear(mpg_d_phi); 
	*/

	//go through all bins
	for (int j = 0; j < m; j++)
	{
		//part of odds equation
		om1 += logFacts[n[j]];
	}
	//final parts of odds equation
	om1 += logFacts[m-1]-logFacts[length+m-1]+((double)length)*log(m);
	return om1;
}

//Equation from gregory and loredo paper to calcluate total odds
//ratio
double log_odds_ratio(double *counts, int length, int m_max, 
					  double nu, double nudot, bool verbose)
{
	//normalize the counts with item 0 at t=0.0s
	normalize_counts(counts, length);
	//the following assumes the counts are ordered
	double t_max = counts[length-1];
	//The total odds ratio
	double odds = 0;
	//go through all possible m values
	for (int m = 2; m <= m_max; m++)
	{
		if (verbose)
			printf("Testing %d-binned model\n",m);
		//Add the next om1 value to the total odds ratio.
		//We also have to remove the log

		odds += exp(log_m_odds_ratio(counts,length,m,nu,
									 nudot,t_max));
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
