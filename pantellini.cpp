#include <iostream>
#include <cmath>
//#include "Random64.h"
#include "pantellini.h"
#include <iomanip>
using namespace std;
//
Crandom::Crandom(unsigned long long j){
    v=4101842887655102017LL; w=1;
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }
unsigned long long Crandom::int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }

double Crandom::exponencial(float tau){
  return -tau*log(r());
}
double Crandom:: gauss(float mu,float sigma){
  return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
}
//
void Ball::InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0)
{i=i0; z=z0; Vx=Vx0; Vy=Vy0; Vz=Vz0;}

void Ball::MovementEq(void)
{
  az = -g;
}

void Ball::Integrate(double dt)
{
  z += Vz*dt+az*dt*dt*0.5;
  Vz += az*dt;
}

void Ball::PrintBall(void)
{
  cout<<" , "<<
    0<<"+"<<R0*0.1<<"*cos(t),"<<
    z<<"+"<<R0*0.1<<"*sin(t)"; 
}


Collider::Collider(void)
{
  for(int i=0;i<N;i++)
    t[i] = 0;
}

void Collider::GetCollisionTime(Ball & ball1, Ball & ball2)
{
  double t12=0,
    Vz1=ball1.Vz, Vz2=ball2.Vz,
    z1=ball1.z,z2=ball2.z, V12 = Vz2-Vz1;
  int i = ball2.i;

  if(V12 == 0.0)
    t12 = 1e10;
  else
    {     
      t12 =abs(z2-z1)/abs(V12);    
      if(t12 > 0 && t12 < tmin)
	{
	  tmin = t12;
	  I = i;
	}
    }
  t[i] = t12;
}

void Collider::CollisionTimeBound(Ball & ball0)
{
  double t0=0, z0 = ball0.z, Vz0 = ball0.Vz;
  t0 = (Vz0 + sqrt(Vz0*Vz0+2*g*z0))/g;
  tmin = t0;
  I = 0;
  t[0] = t0;
}

void Collider::CollideBalls(Ball & ball1, Ball & ball2, Crandom ran, double dt)
{
  double phi = 2*M_PI*ran.r(), P = ran.r();
  double theta = acos(sqrt(P));

  double Vx1 = ball1.Vx, Vy1 = ball1.Vy, Vz1 = ball1.Vz;
  double Vx2 = ball2.Vx, Vy2 = ball2.Vy, Vz2 = ball2.Vz;

  double Ux = Vx2-Vx1, Uy = Vy2-Vy1, Uz = Vz2-Vz1;
  double  U = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);

  double Vx1p = (Vx1+Vx2+U*nx)*0.5, Vx2p = (Vx1+Vx2-U*nx)*0.5;
  double Vy1p = (Vy1+Vy2+U*ny)*0.5, Vy2p = (Vy1+Vy2-U*ny)*0.5;
  double Vz1p = (Vz1+Vz2+U*nz)*0.5, Vz2p = (Vz1+Vz2-U*nz)*0.5;

  ball1.Vx=Vx1p; ball1.Vy=Vy1p; ball1.Vz=Vz1p;
  ball2.Vx=Vx2p; ball2.Vy=Vy2p; ball2.Vz=Vz2p;
  //Integrar...
  //ball1.z+=ball1.Vz*dt; ball2.z+=ball2.Vz*dt;
}
void Collider::CollideBound(Ball & ball0)
{
  double Vz0 = (-1)*ball0.Vz;
  ball0.Vz=Vz0;
}

void Collider::InitPositions(Ball * balls, Crandom ran)
{
  double zi=R0;
  double vix = 0, viy = 0, viz = 0;
  for(int i = 0; i<N; i++)
    {
      zi += R0;
      //viz = ran.r();
      balls[i].InitBall(i, zi, vix, viy, viz);
    }
}

void Collider::MeasureDensity(Ball * balls)
{
  int j = 0;
  double dz = 1e-3, zi = 0, ball_zi;
  int M = h / dz;
  int * n = new int[M];
  int sum_n = 0;   
  for(int l = 0;l<M;l++)
    {
      zi = l*dz;
      for(int i = 0; i<N; i++)
	{
	  ball_zi = balls[i].z;
	  if(ball_zi >= zi && ball_zi < zi+dz)
	    n[l]++;
	}
      if(n[l] != 0)
	cout << zi << "\t" << n[l] << "\n";
      sum_n += n[l];
      if(sum_n >= N)
	break;
    }
  
  delete n;
}

bool Collider::AreBallsSorted(Ball * balls)
{
  bool sorted = true;
  for(int i=0;i<N-1;i++)
    {
      sorted = balls[i].z <= balls[i+1].z;
      if(!sorted)
	break;
    }
  return sorted;
}

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'pelicula.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set yrange[-1:"<<h*0.25<<"]"<<endl;
  cout<<"set xrange["<<-R0<<":"<<R0<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}

void InicieCuadro(void){
  cout<<"plot "<<-R0<<",0 ";
    cout<<" , "<<R0/7<<"*t,0";        //pared de abajo
    cout<<" , "<<R0/7<<"*t,"<<h;     //pared de arriba
    cout<<" , 0,"<<h/7<<"*t";        //pared de la izquierda
    cout<<" , "<<R0<<","<<h/7<<"*t"; //pared de la derecha
}

void TermineCuadro(void){
    cout<<endl;
}


