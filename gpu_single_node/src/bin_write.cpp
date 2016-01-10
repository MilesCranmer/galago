#include <fstream> //for writing the file
#include <stdio.h> //simple output
//this program writes a binary file list of doubles from an array
//Arguments are a pointer to the list in question, the length of the
//list, and the filename

//for fstream, ios, etc.
using namespace std;

void bin_write(double *list, int length, char filename[])
{
	//Open the file to be written in binary
	fstream file(filename, ios::out | ios::binary);
	//Go through all elements of list
	for (int i = 0; i < length; i ++)
	{
		//write the float to the file as binary
		file.write((char*) &list[i], sizeof(double));
	}
}