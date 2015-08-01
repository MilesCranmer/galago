#include <math.h> //for math
//#include <gmpxx.h> //for precision calculation
#include <vector> //to hold search results
#include <stdio.h>
#include <algorithm> //compute max of vector
#include <numeric> //compute sum of vector (accumulate)
#include <omp.h>
//speed not important for final statistics, so optimising this is silly
#define PI 3.14159265359
//this file contains the main statistics functions
//This file relies on main's importing of the logFacts data file

//this number is the max number of threads in a block on
//guillimin
#define BLOCK_SIZE = 1024

using namespace std;

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
