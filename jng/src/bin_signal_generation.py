import numpy as np
import random as rn
import struct
import histogram

def normalize(counts):
	t_min = counts[0]
	return [x-t_min for x in counts]

def sqr(n): return n*n

#non-integral range
#This is inclusive!!!
def frange(start,end,step=1.0):
	curr = start
	#inclusivity
	while curr <= end:
		yield curr	
		curr+=step

def noise(obsTime, ppks):
	t_curr = 0.
	#counts
	counts = []
	#calculate number of photons
	photons = int(obsTime*ppks/1000)
	#calculate 
	for x in range(photons):
		t_curr += rn.expovariate(ppks/1000.)
		if t_curr < obsTime:
			counts.append(t_curr)
	return counts

#This should (eventually) include all the parameters of a neutron 
#star using pulsar catalogues to find information.
#period= the period of the pulsar(s)
#obsTime= total time of the observation(s)
#ppks= photons per kilosecond(ks^-1)
#stdev= standard deviation of each peak(s)
def lgm(period, obsTime, ppks, stdev):
	#Currently these are more artificial variables!!!
	#This will use random.gauss(mu, sigma), which will
	#generate random numbers according to 
	
	#photon counts
	counts = []
	#photons per period (number of photons to generate each loop)
	ppp=ppks/1000.*period
	#Let the first peak occur halfway through the
	#first period.
	for mu in frange(period/2,obsTime,period):
		#Find out how many photons are created in this period:

		#Get number of photons using poisson distribution:
		photons=np.random.poisson(ppp)

		#Create that many photons. Don't record if not in observation.
		for x in range(photons):
			t=np.random.normal(mu,stdev)
			if t>0 and t<obsTime:
				counts.append(t)

	return counts


counts = lgm(10, 1000000, 10, 1)
counts = normalize(sorted(counts))

histogram.dispHisto(counts, 100)
#name of the file
bin_filename = "../data/simple_signal.bin"

#write to file as binary
bin_file = open(bin_filename, "wb")

#iterate through the data
for x in counts:
	#convert to writable binary
	bin = struct.pack('d',x)
	#write this float to the file
	bin_file.write(bin)

#close file
bin_file.close()