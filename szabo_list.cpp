#include <random>
#include <iostream>
#include <fstream>
#include <cmath>

#define N_particles 49

using namespace std;

int cell_index(double ix, double iy){
	//calculate cell index
}

void map_cells(){
	//to list out cell neighbours - Called at the beginning
	
}

void links(){
	//set up head-of-chain array and linked list - Called at each timestep to set up linked list
}

void calc_force(){
	//calculate inter-particle force usig linked list of cells - Called immediately after linked list is set up
}

int main(int argc, char const *argv[]){
	double* head_list, cell_list;
	double Req_mean, R0_mean, L; //INSERT VALUES HERE
	Req_mean = 5.0/6; R0_mean = 1;
	L = Req_mean * sqrt(pi*((double)(N_particles))/rho);

	const int M = (int)(L/R0_mean);
	const int N_cell = M*M*M;
	head_list = new double[M];
	cell_list = new double[N_particles];

	return 0;
}