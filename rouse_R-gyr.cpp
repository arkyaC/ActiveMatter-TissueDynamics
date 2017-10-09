//MAKE SURE THAT THERE IS A data/Rouse/R_gyr FOLDER inside the current directory for storing data
//-----------plotting radius of gyration as a function of monomer units-----------

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <sstream>

#define N_steps 5000
#define default_ensemble_size 200

using namespace std;

random_device rd;
mt19937 gen(rd());//seeding
normal_distribution<> dis(0.0,1.0);//Gaussian random number

void solver(int ctr,int beads) {
  const int N_beads = beads;
  //const int N_beads=5;
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
  
  for(int k = 0;k<N_steps;k++){    
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

  ofstream dump_gyr;
  fileName = "./data/Rouse/R_gyr/gyr_dump_"+to_string(ctr)+".txt";
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

int ensemble(int beads){
  system("exec rm -rf ./data/Rouse/R_gyr/*");//emptying the folder containing old data files

  //float progress = 0.0;
  //int progBarWidth = 60;

  for(int ctr=0;ctr<default_ensemble_size;ctr++){
    /*
    //progress bar
    progress = ((ctr+1)*1.0)/default_ensemble_size;
    int pos = progBarWidth * progress;
    cout<<"\r"<<(progress*100)<<"% complete: ";
    cout<<string(pos,'|');
    cout.flush();
    */
    solver(ctr,beads);
  }
  //cout<<endl;
}

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

int main()
{
  int min_beads = 5,max_beads = 50,N_entries;
  N_entries = max_beads - min_beads + 1;
  double** R_gyr_entries = new double*[N_entries];
  for(int i=0;i<N_entries;i++){
    R_gyr_entries[i] = new double[2];
    R_gyr_entries[i][0] = i + min_beads;
  }
  string fileName, step_data;
  int count;
  ifstream infile;

  for(int ctr_beads = min_beads;ctr_beads<=max_beads;ctr_beads++){

    ensemble(ctr_beads); //calling the function that runs simulation over an ensemble of size default_ensemble_size
    //the folder ./data/Rouse/R_gyr/ now contains gyr_dump*.txt files for no.of beads=ctr_beads

    if(ctr_beads==min_beads){
      //counting how many time steps are covered in each file
      count = 0;
      fileName = "./data/Rouse/R_gyr/gyr_dump_0.txt";
      infile.open(fileName);
      if (! infile){
        cout << "Cannot open input file 0.\n";
        return 1;
      }
      while(getline(infile,step_data)){
        count++;
      }
      infile.close();infile.clear();
    }

    double** avgR_gyr;
    avgR_gyr = new double*[count];
    for(int i=0;i<count;i++){
        avgR_gyr[i] = new double[2];
        avgR_gyr[i][1] = 0.0;
      }
    //evaluating <R_gyr>
    for(int i=0;i<default_ensemble_size;i++){
      fileName = "./data/Rouse/R_gyr/gyr_dump_"+to_string(i)+".txt";
      infile.open(fileName);
      for(int j=count-500;j<count;j++){ //calculating <R_gyr> only for last 500 time steps
        string line_data;
        getline(infile,line_data);
        if(i==0)
          avgR_gyr[j][0] = atof(picker(line_data,1).c_str());
        avgR_gyr[j][1] += (1.0/default_ensemble_size) * atof(picker(line_data,2).c_str());
      }
      infile.close();infile.clear();
    }
    //calculate steady state Radius of gyration <<R_gyr>>
    double steady_state_R = 0;
    for(int i=1;i<=500;i++){ //averaging over last 500 steps in the steady state
      steady_state_R += (1.0/500)*avgR_gyr[count-i][1];
    }
    R_gyr_entries[ctr_beads - min_beads][0] = ctr_beads;
    R_gyr_entries[ctr_beads - min_beads][1] = steady_state_R;

    for(int i=0;i<count;i++){
      delete avgR_gyr[i];
    }
    delete avgR_gyr;
  }

  ofstream dump_gyrVmon;
  fileName = "./data/Rouse/R_gyr-mon.txt";
  dump_gyrVmon.open(fileName);
  for (int i=0;i<N_entries;i++){
    dump_gyrVmon<<R_gyr_entries[i][0]<<"\t"<<R_gyr_entries[i][1]<<"\n";
  }
  dump_gyrVmon.close();

  for(int i=0;i<count;i++){
    delete R_gyr_entries[i];
  }
  delete R_gyr_entries;
}
