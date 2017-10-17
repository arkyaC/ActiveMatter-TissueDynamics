//------------------don't modify delT carelessly, don't lower it from current value-----------------
// SQUARE of end to end distance
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#define N_steps 10000
#define default_ensemble_size 200

using namespace std;

void solver(int ctr) {
  const int N_beads = 10;
  string fileName;

  int N_eff = N_beads + 2;
  double L = 1.0; //length of chain
  double D = 1.0, k_eff = 1.0;//D=diffusion coefficient; define k_eff=k/zeta, zeta being the drag coefficient
  double delT = .1;

  double* xpos = new double[N_eff];
  double* ypos = new double[N_eff];

  random_device rd;
  mt19937 gen(rd());//seeding
  normal_distribution<> dis(0.0,1.0);//Gaussian random number

  for (int i=1;i<=N_beads;i++){
    xpos[i] = (L/(N_beads-1))*(i-1); //initially on a 1D lattice
    ypos[i] = 0;
  }
  xpos[0] = xpos[1]; ypos[0] = ypos[1];//"ghost" nodes
  xpos[N_beads+1] = xpos[N_beads]; ypos[N_beads+1] = ypos[N_beads];

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

  double* e2eDist = new double[N_steps+1];

  for(int k = 0;k<N_steps;k++){

    for(int i=1;i<=N_beads;i++){
      double rhs_i = (-1)*(k_eff)*(2*solX[k][i] - solX[k][i+1] - solX[k][i-1]) ;
      solX[k+1][i] = solX[k][i] + rhs_i * delT + dis(gen)*sqrt(4*D*delT);

      rhs_i = (-1)*(k_eff)*(2*solY[k][i] - solY[k][i+1] - solY[k][i-1]);
      solY[k+1][i] = solY[k][i] + rhs_i * delT + dis(gen)*sqrt(4*D*delT);
    }
    
    e2eDist[k] = pow(solX[k][N_beads]-solX[k][1],2) + pow(solY[k][N_beads]-solY[k][1],2);

    //boundary conditions
    solX[k+1][0] = solX[k+1][1]; solX[k+1][N_beads+1] = solX[k+1][N_beads];
    solY[k+1][0] = solY[k+1][1]; solY[k+1][N_beads+1] = solY[k+1][N_beads];
  }

  ofstream dump_e2e;
  fileName = "./data/Rouse/com_dump_"+to_string(ctr)+".txt";
  dump_e2e.open(fileName);
  for (int i=0;i<N_steps;i++){
    if (i%1==0){
      dump_e2e<<(i+1)<<"\t"<<e2eDist[i]<<"\n";
    }
  }
  dump_e2e.close();

  delete xpos;delete ypos;delete e2eDist;
  for(int i=0;i<N_steps+1;i++){
    delete solX[i];
    delete solY[i];
  }
  delete solX;delete solY;
}

int main(int argc, char const *argv[]){
  int ensmblMax;
  if(argc!=2){
    ensmblMax = default_ensemble_size;
  }
  else{
    ensmblMax=atoi(argv[1]);
  }
  if(ensmblMax<=0)
    ensmblMax = default_ensemble_size;

  system("exec rm ./data/Rouse/*");//emptying the folder containing old data files

  float progress = 0.0;
  int progBarWidth = 60;
  for(int ctr=0;ctr<ensmblMax;ctr++){
    //progress bar
    progress = ((ctr+1)*1.0)/ensmblMax;
    int pos = progBarWidth * progress;
    cout<<"\r"<<(progress*100)<<"% complete: ";
    cout<<string(pos,'|');
    cout.flush();

    solver(ctr);
  }
  cout<<endl;

  return 0;
}