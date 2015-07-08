from numpy import log
import struct #converting to writable binary format

#file to write
bin_filename = "../data/log_facs_2.bin"

#write to file as binary
bin_file = open(bin_filename, "wb")

#print some info about the fits to the terminal
curr = 0

bin = struct.pack('d',curr)
#write this float to the file
#two zeros to start off
bin_file.write(bin)
bin_file.write(bin)

#now go up to a max iterator(around desired number of photons)
for x in range(2, 1000000):
	#next log fact
	curr += log(float(x))
	#convert to writable binary
	bin = struct.pack('d',curr)
	#write this float to the file
	bin_file.write(bin)

	if x%10000==0:
		print x

#close file
bin_file.close()