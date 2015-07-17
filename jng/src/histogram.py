import matplotlib.pyplot as plt

#Display a histogram of the data, with the specified number of bins
def dispHisto(data, bins):
	#Create histogram with specified bins
	plt.hist(data,bins,range=(0,max(data)),color='g',log=False)
	#Other plot preferences
	plt.xlabel('time (s)')
	plt.ylabel('counts')
	plt.title('Histogram of photon counts')
	plt.grid(True)
	plt.show()


#Saves a histogram of the data, with the specified number of bins
def saveHisto(data, bins,filename,results):
	#Create histogram with specified bins
	plt.hist(data,bins,range=(0,max(data)),color='g',log=False)
	#Other plot preferences
	plt.xlabel('time (s)')
	plt.ylabel('counts')
	plt.title('Histogram of photon counts')
	plt.grid(True)
	fig = plt.gcf()
	plt.show()
	plt.draw()
	fig.savefig(filename,dpi=100)