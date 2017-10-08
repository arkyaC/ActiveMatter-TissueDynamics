//MAKE SURE THAT THERE IS A data/Rouse/ FOLDER inside the current directory for storing data
//-----------plotting radius of gyration as a function of monomer units-----------


//------------------------!!!!!!INCOMPLETE!!!!!!!--------------------------------


#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#define N_steps 10000
#define default_ensemble_size 100

using namespace std;

random_device rd;
mt19937 gen(rd());//seeding
normal_distribution<> dis(0.0,1.0);//Gaussian random number

void solver(int ctr) {
  //const int N_beads = beads;
  const int N_beads=5;
  string fileName;

  int N_eff = N_beads + 2;
  double L = 1.0; //length of chain
  double D = 1.0, k_eff = 1.0;//D=diffusion coefficient; define k_eff=k/zeta, zeta being the drag coefficient
  double delT = .1;

  double* xpos = new double[N_eff];
  double* ypos = new double[N_eff];

  /*
  random_device rd;
  mt19937 gen(rd());//seeding
  normal_distribution<> dis(0.0,1.0);//Gaussian random number
  */

  double* comX = new double[N_steps+1];
  double* comY = new double[N_steps+1];
  comX[0]=0.0;comY[0]=0.0; //calculating center of mass position

  for (int i=1;i<=N_beads;i++){
    xpos[i] = (L/(N_beads-1))*(i-1); //initially on a 1D lattice
    comX[0] += xpos[i];
    ypos[i] = 0;
  }
  comX[0] /= N_beads;
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

  double* delcomX = new double[N_steps+1];
  double* delcomY = new double[N_steps+1];
  double* R_gyr = new double[N_steps+1];

  //float progress = 0.0;
  //int progBarWidth = 60;
  
  for(int k = 0;k<N_steps;k++){
    /*
    //progress bar
    if(k%100==0){
      progress = (k*1.0)/N_steps;
      int pos = progBarWidth * progress;
      cout<<"\r"<<(progress*100)<<"% complete: ";
      cout<<string(pos,'|');
      cout.flush();
    }
    */
    
    delcomX[k]=0;delcomY[k]=0;

    for(int i=1;i<=N_beads;i++){
      double rhs_i = (-1)*(k_eff)*(2*solX[k][i] - solX[k][i+1] - solX[k][i-1]) ;
      solX[k+1][i] = solX[k][i] + rhs_i * delT + dis(gen)*sqrt(4*D*delT);
      delcomX[k] += solX[k+1][i] - solX[k][i];

      rhs_i = (-1)*(k_eff)*(2*solY[k][i] - solY[k][i+1] - solY[k][i-1]);
      solY[k+1][i] = solY[k][i] + rhs_i * delT + dis(gen)*sqrt(4*D*delT);
      delcomY[k] += solY[k+1][i] - solY[k][i];
    }
    comX[k+1] = comX[k] + delcomX[k]/N_beads;
    comY[k+1] = comY[k] + delcomY[k]/N_beads;

    R_gyr[k]=0;//initialisation
    for(int i=1;i<=N_beads;i++){
      R_gyr[k] += (1.0/N_beads)*(pow(solX[k][i]-comX[k],2) + pow(solY[k][i]-comY[k],2));
    }
    R_gyr[k] = sqrt(R_gyr[k]);
    
    //boundary conditions
    solX[k+1][0] = solX[k+1][1]; solX[k+1][N_beads+1] = solX[k+1][N_beads];
    solY[k+1][0] = solY[k+1][1]; solY[k+1][N_beads+1] = solY[k+1][N_beads];
  }
  //cout<<endl;//for the progress bar

  ofstream dump_gyr;
  fileName = "./data/Rouse/gyr_dump_"+to_string(ctr)+".txt";
  dump_gyr.open(fileName);
  for (int i=0;i<N_steps;i++){
    if (i%10==0){
      dump_gyr<<(i+1)<<"\t"<<R_gyr[i]<<"\n";
    }
  }
  dump_gyr.close();

  delete xpos;delete ypos;delete comX;delete comY;delete R_gyr;delete delcomX;delete delcomY;
  for(int i=0;i<N_steps+1;i++){
    delete solX[i];
    delete solY[i];
  }
  delete solX;delete solY;
}
//MAIN function
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

  system("exec rm -rf ./data/Rouse/*");//emptying the folder containing old data files

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
