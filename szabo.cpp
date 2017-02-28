#include <random>
#include <iostream>
#include <fstream>
#define pi 3.14159
#define N_steps 5000
#define N_particles 1000
using namespace std;
//producing NaN's (for order parameter) after 800 steps or so
int main(int argc, char const *argv[]) {

  float mu = 1, tau = 1, Req = 5/6, R0 = 1, Fadh = 0.75/6, Frep = 30/6, v0 = 1;
  float noise = 0.6, rhoNorm = 0.3; //tunable parameters
  float L = sqrt(N_particles/(2*rhoNorm)); //rhoMax = 2

  random_device rd;  //Will be used to obtain a seed for the random number engine
  mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  uniform_real_distribution<> dis(0,1);

  float* xpos = new float[N_particles];
  float* ypos = new float[N_particles];
  float* theta = new float[N_particles];
  //cout<<"Test 1"<<endl;
  for (int i=0;i<N_particles;i++){
    xpos[i] = L * dis(gen);
    ypos[i] = L * dis(gen);
    theta[i] = 2 * pi * dis(gen) - pi;
    //cout<<xpos[i]<<"\t"<<ypos[i]<<endl;
  }
  //cout<<"Test 2"<<endl;

  //int N_steps = 1000;
  float delT = 0.05*R0/v0;
  float noiseAmp = noise/sqrt(delT);
  float** solX = new float*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solX[i] = new float[N_particles];
  }
  float** solY = new float*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solY[i] = new float[N_particles];
  }
  //cout<<"Test 2.25"<<endl;
  float** solTheta = new float*[N_steps+1];//[N_particles];
  for(int i=0;i<N_steps+1;i++){
    solTheta[i] = new float[N_particles];
  }
  //cout<<"Test 2.75"<<endl;
  float* order = new float[N_steps];//order parameter of the kth step
  for (int i = 0;i<N_particles;i++){
    //cout<<i<<endl;
    solX[0][i] = xpos[i];
    solY[0][i] = ypos[i];
    solTheta[0][i] = theta[i];
  }
  //cout<<"Test 3"<<endl;

  float* rhsX = new float[N_particles];
  float* rhsY = new float[N_particles];
  float* rhsTheta = new float[N_particles];
//main solver loop
  for(int k = 0;k<N_steps;k++){
    for (int i = 0;i<N_particles;i++){
      float interaxn[2] = {0,0};
      for (int j = 0;j<N_particles;j++){
        if (i==j)
          continue;
        float dx = solX[k][i] - solX[k][j];
        float dy = solY[k][i] - solY[k][j];
        if (abs(dx) > 0.5*L)
          dx = dx - copysign(L,dx);
        if (abs(dy) > 0.5*L)
          dy = dy - copysign(L,dy);
        float dijSq = dx*dx + dy*dy;
        float dij = sqrt(dijSq);
        float eij[2] = {dx/dij,dy/dij};
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
    float avgvel[2] = {0,0};
    for (int i = 0;i<N_particles;i++){
      float ni[2] = {cos(solTheta[k][i]),sin(solTheta[k][i])};
      float vi[2] = {rhsX[i],rhsY[i]};
      float normvi = sqrt(pow(vi[0],2)+pow(vi[1],2));
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
