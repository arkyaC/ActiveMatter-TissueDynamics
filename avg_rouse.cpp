//-----------------to compute ensemble avg----------------
//------------------------INCOMPLETE----------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <sstream>

#define N_files 500
using namespace std;

string picker (string a, int b) { //In today,is,a,good,day if b=3, a will be picked (if delimiter was ,)
  istringstream ss(a);
  string pick;
  int k = 1;
  while(getline(ss, pick, '\t') && k <= b) { //delimiter is \t
      if (k == b){
        return pick;
      }
      k = k + 1;
  }
}

int main(){
	string fileName;
	ifstream infile;
	int count = 0;
	double** avgRcm;
	string step_data;
	
	fileName = "./data/Rouse/com_dump_0.txt";
	infile.open(fileName);
	if (! infile){
		cout << "Cannot open input file.\n";
		return 1;
	}
	while(getline(infile,step_data)){
		count++;
	}
	infile.close();infile.clear();

	avgRcm = new double*[count];
	for(int i=0;i<count;i++){
    	avgRcm[i] = new double[2];
    	avgRcm[i][1] = 0.0;
  	}

	for(int i=0;i<N_files;i++){
		fileName = "./data/Rouse/com_dump_"+to_string(i)+".txt";
		infile.open(fileName);
		for(int j=0;j<count;j++){
			string line_data;
			getline(infile,line_data);
			if(i==0)
				avgRcm[j][0] = atof(picker(line_data,1).c_str());
			avgRcm[j][1] += (1.0/N_files) * atof(picker(line_data,2).c_str());
		}
		infile.close();infile.clear();
	}

	ofstream dump_avg;
	fileName = "./data/Rouse/avg_com.txt";
	dump_avg.open(fileName);
	for (int i=0;i<count;i++){
		dump_avg<<avgRcm[i][0]<<"\t"<<avgRcm[i][1]<<"\n";
	}
	dump_avg.close();	
}