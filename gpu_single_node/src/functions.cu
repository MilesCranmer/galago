#include <thrust/advance.h>
#include <thrust/system_error.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <math.h> //for math
//#include <gmpxx.h> //for precision calculation
#include <vector> //to hold search results
#include <stdio.h>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <algorithm> //compute max of vector
#include <numeric> //compute sum of vector (accumulate)
#include <omp.h>
#include "structures.h"
//speed not important for final statistics, so optimising this is silly
#define PI 3.14159265359
//this file contains the main statistics functions
//This file relies on main's importing of the logFacts data file

//this number is the max number of threads in a block on
//guillimin
#define BLOCK_SIZE = 1024

//extern tells the compiler that this variable is global and
//is already initialized elsewhere
extern double *logFacts;
extern int maxFact;

using namespace std;

//Forward declaration for use in class
double log_odds_ratio(double *counts, int length, int m_max, 
					  double nu, double nudot, bool verbose);


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

/*
//CUDA kernel to create n_mvals*length matrix of bins
__global__ void create_binnings(double *counts, int *mvals,
								int n_mvals, double nu, double nudot,
								unsigned char **binning)
{
	//threads=length
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	double t = counts[idx];
	unsigned char tmp_bin = 0;
	for (int i = 0; i < n_mvals; i++)
	{
		tmp_bin = (unsigned char)((int)(fmod(t*(nu+0.5*t*nudot),1)*mvals[i]));
		binning[i][idx] = tmp_bin;
		binning[i][idx] = 54;

	}
	//n[(int)(fmod(counts[i]*(nu+0.5*counts[i]*nudot),1)*m)]++;
}
*/
//function gets decimal portion of double
__device__ double get_decimal (double x) {return x - (int)x;}

__global__ void create_binnings(double *counts, int *mvals,
								int length,
								int n_mvals, double nu, double nudot,
								unsigned char *binning)
{
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx < length)
	{
		double t = counts[idx];
		int index = idx;
		unsigned char tmp_bin = 0;
		for (int i = 0; i < n_mvals; i++)
		{
			tmp_bin = (int)(get_decimal(t*(nu+0.5*t*nudot))*mvals[i]);
			binning[index] = tmp_bin;
			index += length;
		}
	}
}

//function makes CUDA calls
unsigned char *get_bins(double *counts_d, int length, double *counts_h,
						  int *mvals_d, int *mvals_h, int n_mvals, double nu, double nudot)
{
	unsigned char *binning_h;
	unsigned char *binning_d;
	//initialize thrust arrays
	binning_h = new unsigned char [n_mvals*length];
	binning_h[0] = 100;
	cudaError_t error;
	cudaMalloc((void**)&counts_d,length*sizeof(double));		
	cudaMalloc((void**)&mvals_d,n_mvals*sizeof(double));
	error = cudaMalloc((void**)&binning_d,n_mvals*length);
	if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
	printf("Copying data...\n");
	cudaMemcpy(counts_d,counts_h,length*sizeof(double),
	 		   cudaMemcpyHostToDevice);
	cudaMemcpy(mvals_d,mvals_h,n_mvals*sizeof(double),
	 		   cudaMemcpyHostToDevice);
	error = cudaGetLastError();
	if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
	else {printf("Memory Allocated!\n");}
	printf("Binning data...\n");
	create_binnings<<<40285,1024>>>(counts_d, mvals_d, length, n_mvals, nu, nudot, binning_d);
	thrust::sort(counts_d,counts_d+length);
	//error = cudaThreadSynchronize();	
	if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
	error = cudaMemcpy(binning_h,binning_d,n_mvals*length*sizeof(unsigned char),
			   cudaMemcpyDeviceToHost);
	//error = cudaGetLastError();
	if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
	printf("Done GPU. Cleaning up...\n");
	cudaFree(binning_d);
	cudaFree(counts_d);
	cudaFree(mvals_d);
	return binning_h;
}

/*__global__ void t_bin_counts(thrust::device_vector<double> counts,
							 thrust::device_vector<unsigned char> t_binning,
							 double nu, double nudot,
							 thrust::device_vector<int> mvals)
							 */
__global__ void t_bin_counts(double* counts, int length,
							 unsigned char* t_binning,
							 double nu, double nudot,
							 int* mvals, int n_mvals)
{
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx < length)
	{
		double t = counts[idx];
		int index = idx;
		unsigned char tmp_bin = 0;
		for (int i = 0; i < n_mvals; i++)
		{
			tmp_bin = (unsigned char)(get_decimal(t*(nu+0.5*t*nudot))*mvals[i]);
			t_binning[index] = tmp_bin;
			index += length;
		}
	}
}

__global__ void t_bin_counts_two(double* counts, int length,
							 unsigned char* t_binning,
							 double nu, double nudot)
{
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx < length)
	{
		t_binning[idx] = (unsigned char)(get_decimal(counts[idx]*(nu+0.5*counts[idx]*nudot))*256);
	}
}

__global__ void count_bins(unsigned char *bins, int *histogram_ss, int length)
{

	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx < length-1)
	{
		if (bins[idx] != bins[idx+1])
		{
			//log end of one
			histogram_ss[256+bins[idx]] = idx;
			//log start of other
			histogram_ss[bins[idx+1]] = idx+1;
		}
	}
}
__global__ void get_histo (int *histogram, int *histogram_ss)
{
	int idx = threadIdx.x;
	if (histogram_ss[256+idx] != -1 && histogram_ss[idx] != -1)
		histogram[idx] = histogram_ss[256+idx] - histogram[idx];
}
//function reduces bins by a factor of two
__global__ void reduce_bins_two(int* bins)
{
	int idx = threadIdx.x;	
	bins[idx] = bins[2*idx]+bins[2*idx+1];
}

__global__ void fake_bins(unsigned char *t_binning, int length)
{
	int idx = threadIdx.x;
	t_binning[length+idx] = (unsigned char) idx;
}

__device__ void reduce_bins(int* bins, int m)
{
    for (int i = 0; i < m; i ++)
    {
        bins[i] = bins[2*i]+bins[2*i+1];
    }
}

__global__ void best_five(double *counts, int length,// double *odds_d, 
                          unsigned long long per, double nu_min,
                          double d_nu, double *logFacts_d)
                         // double *nus_d)
{
    //get ID of this core.
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
    //printf("idx=%d\n",idx);
    //make sure last piece
	if (idx < length-1)
	{
        double odds, om1, nudot;
        int m, bins[256]={0};
        //double best[5][3]={0};
        double start = per*d_nu*idx+nu_min;
        double end = start + per*d_nu;
		for (double
		     nu =  start;
			 nu <= end;
			 nu += d_nu)
        {
            odds = 0;
            om1 = 0;
            nudot = 0;
            //nudot=-Pdot/P^2=-v^2*Pdot
            //d_nudot=-nu^2*d_nudot
            //dPdot=Pmin/T^2*P=1/(numax*T^2*nu)
            //for (double
            //        nudot =  nudot_min;
            //        nudot <= nudot_max;
            //        nudot += d_nudot)
            {
                for (int i = 0; i < 256; i++)
                {
                    bins[i] = 0;
                }
                for (int i = 0; i < length; i ++)
                {
                    //With nudot
                    //bins[(unsigned char)(get_decimal(counts[i]*(nu+0.5*counts[i]*nudot))*256)]++;
                    //without nudot
                    bins[(unsigned char)(get_decimal(counts[i]*nu)*256)]++;
                }
                m = 256;
                odds = 0;
                om1 = 0;
                for (int j = 0; j < 256; j++)
                {
                    om1+=logFacts_d[bins[j]];
                }
                om1  += logFacts_d[255]-logFacts_d[length+255]+((double)length)*__log(256.0);
                odds += __exp(om1);
                for (int k = 1; k < 8; k++)
                {
                    m = m >> 1;
                    //make the pointers
                    reduce_bins(bins,m);
                    //	histogram.resize(m);
                    //	binned.resize(m);
                    om1 = 0;
                    //for (int j = 0; j < m; j++)
                    //printf("%d,",binned[j]);
                    //printf("\n");
                    for (int j = 0; j < m; j++)
                    {
                        om1+=logFacts_d[bins[j]];
                    }
                    om1  += logFacts_d[m-1]-logFacts_d[length+m-1]+((double)length)*log((double)m);
                    odds += exp(om1);
                }
                //if (odds > 1e-3)
                odds /= 8;
                odds *= d_nu/nu;
                //put in new best.
				if (odds > 0.1)
				{
                    printf("nu = %e, odds= %e\n", nu, odds);
                    /*
                    //printf("OREONO\n");
					for (int i = 3; i >= 0; i --)
					{
						if (odds < best[i][0])
						{
							for (int j = 3; j >= i + 1; j--) 
							{
								best[j+1][0] = best[j][0];
								best[j+1][1] = best[j][1];
								best[j+1][2] = best[j][2];
							}
							best[i+1][0] = odds;
							best[i+1][1] = nu;
							best[i+1][2] = nudot;
							break;
						}
						else if (i == 0)
						{
							for (int j = 3; j >= 0; j--) 
							{
								best[j+1][0] = best[j][0];
								best[j+1][1] = best[j][1];
								best[j+1][2] = best[j][2];
							}
							best[0][0] = odds;
							best[0][1] = nu;
							best[0][2] = nudot;
						}
					}
                    */
				}
            }
        }
        /*
        for (int i = 0; i < 5; i++)
        {
            odds_d[i+idx*5] = best[i][0];
            nus_d[i+idx*5]  = best[i][1];
        }
        */
    }
    return;
}


double t_odds_two(double *counts_h, int length,
                  double nu_min, double nu_max,
                  double nudot_min, double nudot_max,
                  int verbosity, const char* filename)
{
    try
    {
        //GTX 970 has 1664 cuda cores. 832 per block
        int cores = 1664;
        double d_nu = 1/counts_h[length-1];
        printf("Loading counts...\n");
        thrust::device_vector<double> counts_d(counts_h, counts_h+length);
        double *counts_d_pointer = thrust::raw_pointer_cast(counts_d.data());
        printf("Loading factorials...\n");
        thrust::device_vector<double> logFacts_d(logFacts, logFacts+maxFact);
        double *logFacts_d_pointer = thrust::raw_pointer_cast(logFacts_d.data());

        //calculate number of frequencies to iterate
        unsigned long long op = (unsigned long long)(nu_max-nu_min)/(d_nu);
        unsigned long long per = (unsigned long long)(op/1664.);
        //each thread gets 5 odds allocated to it. Should give the best.
        /*
        thrust::device_vector<double> odds_d(5*cores);
        double *odds_pointer = thrust::raw_pointer_cast(odds_d.data());
        thrust::device_vector<double>   nus_d(5*cores);
        double *nus_pointer = thrust::raw_pointer_cast(nus_d.data());
        */
        printf("Starting search!\n");
        best_five<<<2,1024>>>(counts_d_pointer,length,//odds_pointer,
                per,
                nu_min,d_nu,logFacts_d_pointer);//,nus_pointer);
	    cudaError_t error = cudaThreadSynchronize();	
	    if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
        printf("Search complete\n");
        //printf("Search complete. Retrieving results.\n");
        /*
        thrust::host_vector<double> odds_h = odds_d;
        thrust::host_vector<double>  nus_h =  nus_d;
        for (int i = 0; i < nus_h.size(); i++)
        {
            if (odds_h[i] > 0.1)
                printf("odds = %e, nu = %e\n",odds_h[i],nus_h[i]);
        }
        */
        counts_d.clear();
        counts_d.shrink_to_fit();
        //odds_d.clear();
        //odds_d.shrink_to_fit();
        logFacts_d.clear();
        logFacts_d.shrink_to_fit();
        //nus_d.clear();
        //nus_d.shrink_to_fit();
        /*
        odds_h.clear();
        odds_h.shrink_to_fit();
        nus_h.clear();
        nus_h.shrink_to_fit();
        */
    }
    catch(thrust::system_error &e)
    {
        std::cerr << "Error: "<<e.what() << std::endl;
        exit(-1);
    }
    return 0;
}
				 

double t_odds_two_old(double *counts_h, int length,
				      double nu_min, double nu_max,
				      double nudot_min, double nudot_max,
				      int verbosity, const char* filename)
{
	//the entered mvals should be 2^1 up to 2^8
	try
	{
		double d_nu = 1/counts_h[length-1];
		//double d_nu = 1e-5;
		//double d_nudot = 1e-8;
		//printf("Length: %d\n",length);
		thrust::device_vector<double> counts_d(counts_h, counts_h+length);
		thrust::device_vector<unsigned char> t_binning(length,0);
		thrust::host_vector<int> binned(256,0);
		double *counts_d_pointer = thrust::raw_pointer_cast(counts_d.data());
		unsigned char *t_binning_pointer = thrust::raw_pointer_cast(t_binning.data());
		unsigned int blocks = (unsigned int)(length/1024.0 + 1.0);
		unsigned long counteri = 0;
		double odds = 0;
		double om1 = 0;
		int m;
		int counter = 0;
		double best[5][3] = {0};
		unsigned long opct = (unsigned long)(0.01*(nu_max-nu_min)/d_nu);
		for (double
		     nu =  nu_min;
			 nu <= nu_max;
			 nu += d_nu)
		{
		//	for (double
		//		 nudot =  nudot_min;
		//		 nudot <= nudot_max;
		//		 nudot += d_nudot)
			
			//double d_pdot = d_nu*d_nu/(nu_max*nu);
			double d_nudot = nu*d_nu*d_nu/(nu_max);
			//nudot=-Pdot/P^2=-v^2*Pdot
			//d_nudot=-nu^2*d_nudot
			//dPdot=Pmin/T^2*P=1/(numax*T^2*nu)
			for (double
				 nudot =  nudot_min;
				 nudot <= nudot_max;
				 nudot += d_nudot)
			{

				counteri ++;
				if (counteri >= opct)
				{
					counteri = 0;
					printf("%f percent of the way.\n",100.0*(nu-nu_min)/(nu_max-nu_min));
					ofstream file(filename, ios::app);
					file << "range,"; 
					file << scientific << setprecision(10) << nu_min << "-";
					file << scientific << setprecision(10) << nu << "!";
					file << scientific << setprecision(10) << nu_max << ",";
					file << scientific << nudot_min << "-";
					file << scientific << nudot_max << "\n";
					for (int i = 0; i < 5; i ++)
					{
						//printf("The %dth best odds are %e for a nu of %.9e and nudot -%.9e\n",
						//i+1,best[i][0],best[i][1],best[i][2]);
						file << scientific << best[i][0] << ",";
						file << scientific << setprecision(10) << best[i][1] << ",";
						file << scientific << -best[i][2];
						file << "\n";
					}
					file.close();
				}
				t_bin_counts_two<<<blocks,1024>>>(counts_d_pointer, length, t_binning_pointer, nu, nudot);
				thrust::sort(t_binning.begin(), t_binning.end());
				thrust::device_vector<int> histogram(256,0);

				thrust::counting_iterator<int> search_begin(0);
				thrust::upper_bound(t_binning.begin(), t_binning.end(),
									search_begin, search_begin + 256,
									histogram.begin());
				thrust::adjacent_difference(histogram.begin(), histogram.end(),
										    histogram.begin());
				binned=histogram;
				m = 256;
				odds = 0;
				om1 = 0;
				for (int j = 0; j < 256; j++)
				{
					om1+=logFacts[binned[j]];
				}
				om1  += logFacts[255]-logFacts[length+255]+((double)length)*log(256);
				odds += exp(om1);
				for (int k = 1; k < 8; k++)
				{
					/*
					for (int x = 0; x < m; x ++)
					{
						printf("%d,",binned[x]);
					}
					printf("\n");
					*/
					m = m >> 1;
					//printf("m=%d\n",m);
					//make the pointers
					int *bins_d = thrust::raw_pointer_cast(histogram.data());
					reduce_bins_two<<<1,m>>>(bins_d);
				//	histogram.resize(m);
				//	binned.resize(m);
					binned = histogram;
					om1 = 0;
					//for (int j = 0; j < m; j++)
					//printf("%d,",binned[j]);
					//printf("\n");
					for (int j = 0; j < m; j++)
					{
						om1+=logFacts[binned[j]];
					}
					om1  += logFacts[m-1]-logFacts[length+m-1]+((double)length)*log(m);
					odds += exp(om1);
				}
				//if (odds > 1e-3)
				odds /= 8;
				odds *= d_nu/nu;
				//results.nu.push_back(nu);
				//results.nudot.push_back(nudot);
				//results.odds.push_back(odds);
				//if (counter %50000==0 || odds > 1e-4)
				/*
				if (verbosity == 2 || (verbosity == 1 && odds > 1e-3) || (verbosity == 0 && odds > 1e-1))
				{
					printf("Search %d gives odds of %e for nu %.9e and nudot -%.9e\n",counter,odds,nu,nudot);
				}
				else if (verbosity == 1 && counter%50000==0)
				{
					printf("On search %d, and nu=%.9e Hz\n",counter,nu);	
				}
				*/
				
				if (odds > best[4][0])
				{
					for (int i = 3; i >= 0; i --)
					{
						if (odds < best[i][0])
						{
							for (int j = 3; j >= i + 1; j--) 
							{
								best[j+1][0] = best[j][0];
								best[j+1][1] = best[j][1];
								best[j+1][2] = best[j][2];
							}
							best[i+1][0] = odds;
							best[i+1][1] = nu;
							best[i+1][2] = nudot;
							break;
						}
						else if (i == 0)
						{
							for (int j = 3; j >= 0; j--) 
							{
								best[j+1][0] = best[j][0];
								best[j+1][1] = best[j][1];
								best[j+1][2] = best[j][2];
							}
							best[0][0] = odds;
							best[0][1] = nu;
							best[0][2] = nudot;
						}
					}
				}
			}
		}
		//clear up space
		counts_d.clear();
		counts_d.shrink_to_fit();
		ofstream file(filename, ios::app);
		file << "range,"; 
		file << scientific << setprecision(10) << nu_min << "-";
		file << scientific << setprecision(10) << nu_max << ",";
		file << scientific << nudot_min << "-";
		file << scientific << nudot_max << "\n";
		for (int i = 0; i < 5; i ++)
		{
			//printf("The %dth best odds are %e for a nu of %.9e and nudot -%.9e\n",
				   //i+1,best[i][0],best[i][1],best[i][2]);
			file << scientific << best[i][0] << ",";
			file << scientific << setprecision(10) << best[i][1] << ",";
			file << scientific << -best[i][2];
			file << "\n";
		}
		file.close();
		//results.write_best(10,filename);
		//best[5][3];
		//keep reducing bins
		//int j = results.max_odds_i();
		//printf("\nThe best odds are: %e, which occur for nu of %e Hz and"
			   //" nudot of -%e Hz/s\n\n",
		//		results.odds[j], results.nu[j], results.nudot[j]);
		//printf("%d searches completed\n",counter);
		return 0;
	}
	catch(thrust::system_error &err)
	{
		std::cerr << "Error doing this: " << err.what() << std::endl;
		return 1;
	}
}


double t_odds(double *counts_h, int length,
		double nu, double nudot,
		int *mvals_h, int n_mvals)
{
	try
	{
		thrust::device_vector<double> counts_d(counts_h, counts_h+length);
		thrust::device_vector<unsigned char> t_binning(length*n_mvals,0);
		thrust::device_vector<int> mvals_d(mvals_h, mvals_h+n_mvals);
		double *counts_d_pointer = thrust::raw_pointer_cast(counts_d.data());
		unsigned char *t_binning_pointer = thrust::raw_pointer_cast(t_binning.data());
		int *mvals_d_pointer = thrust::raw_pointer_cast(mvals_d.data());
		t_bin_counts<<<40285,1024>>>(counts_d_pointer, length, t_binning_pointer, nu, nudot, 
									 mvals_d_pointer, n_mvals);
		//clear up space
		counts_d.clear();
		counts_d.shrink_to_fit();
		//iterate through segments of array
		thrust::device_vector<unsigned char>::iterator iter_start = t_binning.begin();
		thrust::device_vector<unsigned char>::iterator iter_end   = t_binning.begin();
		//sort parts of array
		for (int i = 0; i < n_mvals; i++)
		{
			thrust::advance(iter_end,length);
			thrust::sort(iter_start, iter_end);
			//thrust::sort(&t_binning[i*length],&t_binning[i*length + length]);
			thrust::advance(iter_start,length);
		}
		
		double odds = 0;
		iter_start = t_binning.begin();
		iter_end = t_binning.begin();
		for (int i = 0; i < n_mvals; i++)
		{
			thrust::advance(iter_end,length);
			thrust::device_vector<int> histogram(mvals_h[i],0);
			thrust::host_vector<unsigned char> histo_vals_h(mvals_h[i],0);
			for (unsigned char j = 0; j < mvals_h[i]; j++)
			{
				histo_vals_h[j] = j;
			}
			thrust::device_vector<unsigned char> histo_vals=histo_vals_h;
			thrust::reduce_by_key(iter_start, iter_end,	
					thrust::constant_iterator<int>(1),
					histo_vals.begin(),
					histogram.begin());
			thrust::advance(iter_start,length);
			//load these values back to the host, as has been binned
			thrust::host_vector<int> binned = histogram;
			double om1 = 0;
			for (int j = 0; j < binned.size(); j++)
			{
				om1+=logFacts[binned[j]];
			}
			om1  += logFacts[mvals_h[i]-1]-logFacts[length+mvals_h[i]-1]+((double)length)*log(mvals_h[i]);
			odds += exp(om1);
		}
		return odds;
	}
	catch(thrust::system_error &err)
	{
		std::cerr << "Error doing this: " << err.what() << std::endl;
		exit(-1);
	}

}

//function uploads static data to the GPU at start of MPI proc
void upload_data(double *counts_h, double *counts_d, int length,
				 int *mvals_h, int *mvals_d, int n_mvals)
{
	printf("1000 toa = %f\n",counts_h[1000]);
	//cudaMalloc((void**)&counts_d,length*sizeof(double));		
	//cudaMalloc((void**)&mvals_d,n_mvals*sizeof(double));
	//cudaMemcpy(counts_d,counts_h,length*sizeof(double),
	// 		   cudaMemcpyHostToDevice);
	//cudaMemcpy(mvals_d,mvals_h,n_mvals*sizeof(double),
	//		   cudaMemcpyHostToDevice);
	cudaError_t error = cudaGetLastError();
	if (error!=cudaSuccess) {printf("Error! %s\n",cudaGetErrorString(error));}
	else{printf("Static data uploaded!\n");}
}

void free_data(double *counts_d, int *mvals_d)
{
	//cudaFree(mvals_d);
	//cudaFree(counts_d);
}

double bins_to_odds(unsigned char *bins, int length, 
					int *mvals, int n_mvals)
{
	double odds = 0;
	for (int i = 0; i < n_mvals; i ++)
	{
		double om1 = 0;
		int m = mvals[i];	
		int n[m];
		for (int j = 0; j < m; j ++)
		{
			n[j] = 0;
		}
		for (int k = i*length; k < (i+1)*length; k++)
		{
			n[bins[k]]++;
		}
		/*
		#pragma omp parallel
		{
			int ni[m];
			for (int r = 0; r < m; r ++)
			{
				ni[r] = 0;
			}
			#pragma omp for
			for (int k = i*length; k < (i+1)*length; k++)
			{
				ni[bins[k]]++;
			}
			#pragma omp critical
			for (int q = 0; q < m; q++)
			{
				n[q] += ni[q];	
			}
		}
		*/
		for (int l = 0; l < m; l++)
		{
			//part of odds equation
			om1 += logFacts[n[l]];
		}
		om1 += logFacts[m-1]-logFacts[length+m-1]+((double)length)*log(m);
		odds += exp(om1);
	}
	return odds;
}

				

//Equation from gregory and loredo paper to calcluate odds ratio
//of m-binned stepwise model w.r.t. constant model
double log_m_odds_ratio(double *counts, int length, int m, 
					  double nu, double nudot,
					  double t_max)
{
	//create all the bins, init to zero counts
	unsigned int ng[m];

	//init to zero
	for (int j = 0; j < m; j++)
	{
		ng[j] = 0;
	}
	
	//split up into threads
	#pragma omp parallel default(shared)
	{
		//create temp bins for thread
		unsigned int n[m];
		for (int j = 0; j < m; j++)
		{
			n[j] = 0;
		}

		//variables used in binnings
		//gets position in nu
		//long double phi, d_phi;
		//double phi;
		//gets bin
		//int k;
		//bin the photons
		#pragma omp for 
		for (int i = 0; i < length; i++)
		{
			
			//d_phi = 0.5*counts[i]*nudot*counts[i];
			//get position in nu of photon
			//phi = fmod(counts[i]*nu+d_phi,1);
			//get corresponding bin	
			//k = (int)(fmod(counts[i]*(nu+0.5*counts[i]*nudot),1)*m);
			//one more count
			n[(int)(fmod(counts[i]*(nu+0.5*counts[i]*nudot),1)*m)]++;
			
		}

		//combine n values
		#pragma omp critical
		for (int j = 0; j < m; j++)
		{
			ng[j] += n[j];
		}
	}
	//odds to return
	double om1 = 0.0;
	//go through all bins
	for (int j = 0; j < m; j++)
	{
	
		//part of odds equation
		om1 += logFacts[ng[j]];
	}
	//final parts of odds equation
	om1 += logFacts[m-1]-logFacts[length+m-1]+((double)length)*log(m);
	return om1;
}

//Equation from gregory and loredo paper to calcluate total odds
//ratio
double log_odds_ratio(double *counts, int length, int *mvals, int n_mvals, 
					  double nu, double nudot, bool verbose)
{
	//normalize the counts with item 0 at t=0.0s
	normalize_counts(counts, length);
	//the following assumes the counts are ordered
	double t_max = counts[length-1];
	//The total odds ratio
	double odds = 0.0;
	//go through all possible m values
	for (int i = 0; i <= n_mvals; i++)
	{
		if (verbose)
			printf("Testing %d-binned model\n",i);
		//Add the next om1 value to the total odds ratio.
		//We also have to remove the log

		odds += exp(log_m_odds_ratio(counts,length,mvals[i],nu,
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

