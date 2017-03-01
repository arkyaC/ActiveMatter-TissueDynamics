#include <random>
#include <iostream>
#include <fstream>
#include <cmath>
#define pi 3.14159
#define N_steps 10000
#define N_particles 1000
using namespace std;
//producing NaN's (for order parameter) after 800 steps or so
int main(int argc, char const *argv[]) {

  double mu = 1, tau = 1, Req = 5/6, R0 = 1, Fadh = 0.75/6, Frep = 30/6, v0 = 1;
  double noise = 0.6, rhoNorm = 0.3; //tunable parameters
  double L = R0 * sqrt(N_particles/(2*rhoNorm)); //rhoMax = 2

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dis(0,1); //uniform distribution on [0,1)

  double* xpos = new double[N_particles];
  double* ypos = new double[N_particles];
  double* theta = new double[N_particles];
  //cout<<"Test 1"<<endl;
  for (int i=0;i<N_particles;i++){
    xpos[i] = L * dis(gen);
    ypos[i] = L * dis(gen);
    theta[i] = 2 * pi * dis(gen) - pi;
    //cout<<xpos[i]<<"\t"<<ypos[i]<<endl;
  }
  //cout<<"Test 2"<<endl;

  //int N_steps = 1000;
  double delT = 0.05*R0/v0;
  double noiseAmp = noise/sqrt(delT);
  double** solX = new double*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solX[i] = new double[N_particles];
  }
  double** solY = new double*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solY[i] = new double[N_particles];
  }
  //cout<<"Test 2.25"<<endl;
  double** solTheta = new double*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solTheta[i] = new double[N_particles];
  }
  //cout<<"Test 2.75"<<endl;
  double* order = new double[N_steps];//order parameter of the kth step
  for (int i = 0;i<N_particles;i++){
    //cout<<i<<endl;
    solX[0][i] = xpos[i];//Initialize 0th step
    solY[0][i] = ypos[i];
    solTheta[0][i] = theta[i];
  }
  //cout<<"Test 3"<<endl;

  double* rhsX = new double[N_particles];
  double* rhsY = new double[N_particles];
  double* rhsTheta = new double[N_particles];
//main solver loop
  for(int k = 0;k<N_steps;k++){
    for (int i = 0;i<N_particles;i++){
      double interaxn[2] = {0,0}; //interaction force vector for ith particle
      for (int j = 0;j<N_particles;j++){
        if (i==j)
          continue;
        double dx = solX[k][i] - solX[k][j];//dij vector = (dx,dy)
        double dy = solY[k][i] - solY[k][j];
        if (abs(dx) > 0.5*L)
          dx = dx - copysign(L,dx);
        if (abs(dy) > 0.5*L)
          dy = dy - copysign(L,dy);
        double dijSq = dx*dx + dy*dy;
        double dij = sqrt(dijSq);
        double eij[2] = {dx/dij,dy/dij};
        if (dijSq<Req*Req){
            interaxn[0] -= Frep*(dij-Req)/Req * eij[0];
            interaxn[1] -= Frep*(dij-Req)/Req * eij[1];
          }
        else if (dijSq>=Req*Req && dijSq<R0*R0){
            interaxn[0] -= Fadh*(dij-Req)/(R0-Req) * eij[0];
            interaxn[1] -= Fadh*(dij-Req)/(R0-Req) * eij[1];
          }
      }
      interaxn[0] *= mu;
      interaxn[1] *= mu;
      rhsX[i] = v0*cos(solTheta[k][i])+interaxn[0];
      rhsY[i] = v0*sin(solTheta[k][i])+interaxn[1];
    }
    double avgvel[2] = {0,0};
    for (int i = 0;i<N_particles;i++){
      double ni[2] = {cos(solTheta[k][i]),sin(solTheta[k][i])};
      double vi[2] = {rhsX[i],rhsY[i]};
      double normvi = sqrt(pow(vi[0],2)+pow(vi[1],2));
      vi[0] = vi[0]/normvi;
      vi[1] = vi[1]/normvi;
      avgvel[0] += vi[0];
      avgvel[1] += vi[1];
      rhsTheta[i] = asin(ni[0]*vi[1]-ni[1]*vi[0]);
      solX[k+1][i] = solX[k][i] + delT*rhsX[i];
      solX[k+1][i] -= L*floor(solX[k+1][i]/L);//wrapping around excesses
      solY[k+1][i] = solY[k][i] + delT*rhsY[i];
      solY[k+1][i] -= L*floor(solY[k+1][i]/L);//wrapping around excesses
      solTheta[k+1][i] = solTheta[k][i] + delT*(rhsTheta[i]/tau + noiseAmp * (dis(gen)-0.5));
    }
    order[k] = sqrt(pow(avgvel[0],2)+pow(avgvel[1],2))/N_particles;
  }

  ofstream dump_data;
  dump_data.open("dump.txt");
  dump_data<<"Timestep\tOrder Parameter\n";
  for (int i=0;i<N_steps;i++){
    dump_data<<(i+1)<<"\t"<<order[i]<<"\n";
  }
  dump_data.close();
  return 0;
}
