//------------------don't modify delT carelessly-----------------

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <sstream>
#define N_steps 40000
#define pi 3.141593
#define ensemble_size 200

using namespace std;

random_device rd_x;
random_device rd_y;
mt19937 gen_x(rd_x());//seeding
mt19937 gen_y(rd_y());//seeding
normal_distribution<> dis(0.0,1.0);//Gaussian random number
uniform_real_distribution<> dis_uni(-1.0,1.0);

void ensemble(int ctr){
	const int N_beads = 10;

	int N_eff = N_beads + 2;
	double epsilon = 1.0; //strength of WCA; range sigma = b
	double a = 5.0; //equilibrium bond length

	double tau = a/sqrt(epsilon);//WCA time scale
	double kbT = 1.0;//temperature
	double gamma = 1.0;//friction coefficient
	double k_sp = 100.0; //spring constant
	double b = sqrt(2*kbT/k_sp); //Equipartition theorem b^2=<x^2>+<y^2>
	double D = kbT/gamma;

	double delT = .02*gamma/k_sp; 

	//-------------------------setting up the gaussian chain-------------------
	double* xpos = new double[N_eff];
	double* ypos = new double[N_eff];

	xpos[1] = 0.0; ypos[1] = 0.0;
	for (int i=2;i<=N_beads;i++){
		//Gaussian Chain i.c.

			// double r_val = a + b*dis(gen_x);//gaussian in r

			// double theta = pi * dis_uni(gen_y);//uniform in theta
			// // double theta = pi * dis_uni(gen_x);//uniform in theta

			// xpos[i] = xpos[i-1] + r_val*cos(theta); //initially a gaussian chain
			// ypos[i] = ypos[i-1] + r_val*sin(theta); 

		//Linear i.c.
		xpos[i] = (i-1)*a+b*dis(gen_x);
		ypos[i] = 0;
	}
	xpos[0] = xpos[1]; ypos[0] = ypos[1];//"ghost" nodes
	xpos[N_beads+1] = xpos[N_beads]; ypos[N_beads+1] = ypos[N_beads];


	// -------------testing i.c.------------------------
	    // cout<<"xpos array"<<endl;
	    // for(int i=1;i<=N_beads;i++){
	    //   cout<<xpos[i]<<",";
	    // }
	    // cout<<endl;
	    // cout<<"ypos array"<<endl;
	    // for(int i=1;i<N_beads;i++){
	    //   cout<<ypos[i]<<",";
	    // }
	    // cout<<endl;

	//--------------------------------------------------


	//-------------------------rouse dynamics---------------------------------

	//Declare solution matrices
	double** solX = new double*[N_steps+1];
	for(int i=0;i<N_steps+1;i++){
		solX[i] = new double[N_eff];
	}
	double** solY = new double*[N_steps+1];
	for(int i=0;i<N_steps+1;i++){
		solY[i] = new double[N_eff];
	}
	//Initialize 0th step
	for (int i = 0;i<N_eff;i++){
		solX[0][i] = xpos[i];
		solY[0][i] = ypos[i];
	}

	double* e2e_dist = new double[N_steps+1];//stores (end-to-end dist)^2

	//progress bar
		// float progress = 0.0;
		// int progBarWidth = 60;
	for(int k = 0;k<N_steps;k++){
    //progress bar
		// if(k%100==0){
		// 	progress = (k*1.0)/N_steps;
		// 	int pos = progBarWidth * progress;
		// 	cout<<"\r"<<(progress*100)<<"% complete: ";
		// 	cout<<string(pos,'|');
		// 	cout.flush();
		// }

    for(int i=1;i<=N_beads;i++){
      	//spring force for adjacent monomers

		double dist_p = sqrt(pow(solX[k][i] - solX[k][i+1],2)+pow(solY[k][i] - solY[k][i+1],2));
		double dist_m = sqrt(pow(solX[k][i] - solX[k][i-1],2)+pow(solY[k][i] - solY[k][i-1],2));
		double rhs_x=0,rhs_y=0;
		if(i==1){
			rhs_x = -(k_sp/gamma)*(solX[k][i] - solX[k][i+1] - a*(solX[k][i]-solX[k][i+1])/dist_p);
			rhs_y = -(k_sp/gamma)*(solY[k][i] - solY[k][i+1] - a*(solY[k][i]-solY[k][i+1])/dist_p);
		}
		else if(i==N_beads){
			rhs_x = -(k_sp/gamma)*(solX[k][i] - solX[k][i-1] - a*(solX[k][i]-solX[k][i-1])/dist_m);
			rhs_y = -(k_sp/gamma)*(solY[k][i] - solY[k][i-1] - a*(solY[k][i]-solY[k][i-1])/dist_m);
		}
		else{
			rhs_x = -(k_sp/gamma)*(2*solX[k][i] - solX[k][i+1] - solX[k][i-1] - a*(solX[k][i]-solX[k][i+1])/dist_p - a*(solX[k][i]-solX[k][i-1])/dist_m);
			rhs_y = -(k_sp/gamma)*(2*solY[k][i] - solY[k][i+1] - solY[k][i-1] - a*(solY[k][i]-solY[k][i+1])/dist_p - a*(solY[k][i]-solY[k][i-1])/dist_m);
		}

		//WCA potential for non-adjacent monomers
		for(int j=1;j<=N_beads;j++){//need to make this more efficient!!!!
			if(abs(j-i)<2)
			  continue;
			double rijSq = pow(solX[k][i]-solX[k][j],2) + pow(solY[k][i]-solY[k][j],2);
			if(rijSq<pow(2,1.0/3)*a*a){//WCA contribution
			  rhs_x += (1.0/gamma)*24*epsilon*pow(a,6)*(2*pow(a*a/rijSq,3)-1)*(solX[k][i]-solX[k][j])/pow(rijSq,4);
			  rhs_y += (1.0/gamma)*24*epsilon*pow(a,6)*(2*pow(a*a/rijSq,3)-1)*(solY[k][i]-solY[k][j])/pow(rijSq,4);
			}
		}

		solX[k+1][i] = solX[k][i] + rhs_x * delT + dis(gen_x)*sqrt(2*D*delT);

		solY[k+1][i] = solY[k][i] + rhs_y * delT + dis(gen_y)*sqrt(2*D*delT);
		// solY[k+1][i] = solY[k][i] + rhs_y * delT + dis(gen_x)*sqrt(2*D*delT);

    }

    //End-to-End distance square
    e2e_dist[k] = pow(solX[k][N_beads]-solX[k][1],2) + pow(solY[k][N_beads]-solY[k][1],2);

    //boundary conditions
    solX[k+1][0] = solX[k+1][1]; solX[k+1][N_beads+1] = solX[k+1][N_beads];
    solY[k+1][0] = solY[k+1][1]; solY[k+1][N_beads+1] = solY[k+1][N_beads];
  }
  cout<<endl;//for progress bar

  ofstream dump_e2e;
  string fileName = "./data/wca/e2e_dump_"+to_string(ctr)+".txt";
  dump_e2e.open(fileName);
  for (int i=0;i<N_steps;i++){
    // if (i%10==0)
      dump_e2e<<(i+1)<<"\t"<<e2e_dist[i]<<"\n"; //dumping end to end distance in data file
  }
  dump_e2e.close();
}

//-------------------------------This is a string parsing function------------------------------------
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

int main(int argc, char const *argv[]) {
	int ctr = 0, count = 0;

	//progress bar
		float progress = 0.0;
		int progBarWidth = 60;
	for(;ctr<ensemble_size;ctr++){
		//progress bar
			progress = (ctr*1.0)/ensemble_size;
			int pos = progBarWidth * progress;
			cout<<"\r"<<(progress*100)<<"% complete: ";
			cout<<string(pos,'|');
			cout.flush();

		ensemble(ctr);
	}
	string fileName = "./data/wca/e2e_dump_0.txt",step_data;
	ifstream infile;
	infile.open(fileName);
	if (! infile){
    	cout << "Cannot open input file 0.\n";
        return 1;
    }
	while(getline(infile,step_data))
		count++;
	infile.close();infile.clear();

	double* avg_e2e;
    avg_e2e = new double[count];
    for(int i=0;i<ensemble_size;i++){
    	fileName = "./data/wca/e2e_dump_"+to_string(i)+".txt";
    	infile.open(fileName);
    	for(int j=(int)(0.5*count);j<count;j++){
    		string line_data;
	        getline(infile,line_data);
	        avg_e2e[j] += (1.0/ensemble_size) * atof(picker(line_data,2).c_str());
    	}
    	infile.close();infile.clear();
    }

    ofstream dump_e2e;
	fileName = "./data/wca/avg_e2e.txt";
	dump_e2e.open(fileName);
	for (int i=0;i<count;i++){
	// if (i%10==0)
		dump_e2e<<(i+1)<<"\t"<<avg_e2e[i]<<"\n"; //dumping end to end distance in data file
	}
	dump_e2e.close();
	delete avg_e2e;
    
	return 0;
}
