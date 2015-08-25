#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "bin_read.cpp"

using namespace std;

bool is_file_exist(const char *fileName)
{
    ifstream infile(fileName);
    return infile.good();
}

int main()
{
	ifstream file;
	file.open("grids.csv", ios::in);
	string curr_settings;
	if (!file.is_open()) 
	{
		cout << "Error opening file" << endl; 
		return 0;
	}
	//Broadcast max, min of nu, nudot
	//broadcast logfacts
	//
	//send a process a list of times
	//
	//MPI_Procs should now only get
	//
	//
	while (getline(file, curr_settings))
	{
		int x, y;
		x = atoi(curr_settings.substr(0,3).c_str());
		y = atoi(curr_settings.substr(4,3).c_str());
		stringstream ss;
		ss << "/Users/miles/Downloads/11819_package/merged2/grid/grid_";
		ss << x;
		ss << "_";
		ss << y;
		ss << ".bin";
		string filename = ss.str();
		if (is_file_exist(filename.c_str()))
		{
			//cout << filename << endl;
			//char filename[] = "/Users/miles/Downloads/11819_package/merged2/grid/grid_"+to_string(x)+"_"+to_string(y)+".bin";
			int length = bin_size((char*)filename.c_str());
			//int length = bin_size((char*)"/Users/miles/Downloads/11819_package/merged2/grid/grid_313_229.bin");
			printf("File is %d doubles long\n", length);
		}


		//grid_316_227.bin
	}
	file.close();
	return 0;
}