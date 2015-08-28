#include <fstream>//file input
#include <stdio.h>//console output
//this program reads a binary file list of doubles into an array

//for fstream, ios, etc.
using namespace std;

double * bin_read(char filename[])
{
	//first figure out size of file:
	streampos size;
	//open binary file to read
	fstream file(filename, ios::in | ios::binary | ios::ate);
	//make sure file is open
	if (!file.is_open()){printf("Cannot open file\n");return NULL;}
	//Get size of file
	size = file.tellg();
	//declare double of the correct size:
	int doubleLength = size/sizeof(double);
	//need to allocate so memory will persist in main function
	double *list;
	list = new double [doubleLength];
	//Get position at beginning of file
	file.seekg (0, ios::beg);
	//append all doubles to list
	for (int i = 0; i < doubleLength; i++)
	{
		//temporary character which will hold double
		unsigned char tmp[sizeof(double)];
		//reinterpret_cast does not compile to the CPU - compiler 
		//directive to treat sequence of bits in (expr) as if it 
		//had the other type.
		file.read(reinterpret_cast<char*>(tmp), sizeof(double));
		//Add float to array
		if (file.good())
		{
			list[i] = reinterpret_cast<double&>(tmp);
		}

	}
	file.close();
	return list;
}

//This function gets the number of doubles in the file
int bin_size(char filename[])
{
	streampos size;
	//open binary file to read
	fstream file(filename, ios::in | ios::binary | ios::ate);
	//make sure file is open
	if (!file.is_open()){printf("Cannot open file\n");return 0;}
	//Get size of file
	size = file.tellg();
	//declare double of the correct size:
	return size/sizeof(double);
}