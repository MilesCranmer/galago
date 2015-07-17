import astropy.io.fits as fits #fits input
import struct #converting to writable binary format

fits_filename = "../data/B1821.fits"
bin_filename = "../data/B1821_counts.bin"

#write to file as binary
bin_file = open(bin_filename, "wb")

#open hdu
hdulist = fits.open(fits_filename)
#print some info about the fits to the terminal
print hdulist.info()
i = 0
#iterate through the data
for x in hdulist[1].data:
	#load time of arrival
	curr = float(x[0])
	#convert to writable binary
	bin = struct.pack('d',curr)
	#write this float to the file
	bin_file.write(bin)
	i+=1
	#set limit
	if i == 200000:
		break

#close files
hdulist.close()
bin_file.close()