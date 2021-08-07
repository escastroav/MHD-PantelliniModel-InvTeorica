#include <iostream>
#include <cmath>
#include "Random64.h"
#include <iomanip>
#include <cstdlib>
using namespace std;

const double size = 1;
const int N = 50*size;
const double R0 = 1/size;

const double h = N*R0;
const double g = 9.8e-3*R0;
const double H = h*.2; // kT/mg=rmax/5

const double tf = sqrt(H/g);

class Ball
{
private:
  int i = 0;
  double z=0,
    Vx=0,
    Vy=0,
    Vz=0;
  double az=0;
public:
  void InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0);
  double GetZ(void){return z;};void SetZ(double zi){z=zi;};
  double GetVx(void){return Vx;};void SetVx(double Vxi){Vx=Vxi;};
  double GetVy(void){return Vy;};void SetVy(double Vyi){Vy=Vyi;};
  double GetVz(void){return Vz;};void SetVz(double Vzi){Vz=Vzi;};
  int GetI(void){return i;};
  void MovementEq(void);
  void Integrate(double dt);
  void PrintBall(void);
  friend class Collider;
};
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
class Collider
{
private:
  double t[N];
  double tmin=0.0;
  double I = 0;
  double * n = nullptr;
  double * f = nullptr;
  double dz = 0.245*H, dv=0.1;
  int M = 0, Mv = 100;
public:
  Collider(void);
  ~Collider(void);
  double GetTmin(void){return tmin;};
  int GetIndex(void){return I;};
  void GetCollisionTime(Ball & ball1, Ball & ball2);
  void CollisionTimeBound(Ball & ball0);
  void CollideBalls(Ball & ball1, Ball & ball2, Crandom rand, double dt);
  void CollideBallsAlt(Ball & ball1, Ball & ball2, Crandom rand, double dt);
  void CollideBound(Ball & ball0);
  void InitPositions(Ball * balls, Crandom ran);
  void MeasureDensity(Ball * balls);
  void PrintDensity(void);
  double TotalEnergy(Ball * balls);
  bool AreBallsSorted(Ball * balls);
  void VelocityDistribution(Ball * balls);
  void PrintVelocities(void);
};
Collider::Collider(void)
{
  M = 100;
  dz = h / M;
  n = new double[M];
  f = new double[Mv];
  for(int i=0;i<N;i++)
    t[i] = 0;
}
Collider::~Collider(void)
{
  delete n; delete f;
}
void Collider::GetCollisionTime(Ball & ball1, Ball & ball2)
{
  double t12=0,
    Vz1=ball1.Vz, Vz2=ball2.Vz,
    z1=ball1.z,z2=ball2.z, V12 = Vz2-Vz1;
  int i = ball2.i;
  if(V12 >= 0.0)
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
  ran.Reset((unsigned long long)rand());
  double phi = 2*M_PI*ran.r(), P = ran.r();
  double theta = acos(sqrt(P));// cout << theta << "\n";

  double Vx1 = ball1.Vx, Vy1 = ball1.Vy, Vz1 = ball1.Vz;
  double Vx2 = ball2.Vx, Vy2 = ball2.Vy, Vz2 = ball2.Vz;

  double Ux = Vx2-Vx1, Uy = Vy2-Vy1, Uz = Vz2-Vz1;
  double  U = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);

  double Vx1p = (Vx1+Vx2-U*nx)*0.5, Vx2p = (Vx1+Vx2+U*nx)*0.5;
  double Vy1p = (Vy1+Vy2-U*ny)*0.5, Vy2p = (Vy1+Vy2+U*ny)*0.5;
  double Vz1p = (Vz1+Vz2-U*nz)*0.5, Vz2p = (Vz1+Vz2+U*nz)*0.5;

  ball1.Vx=Vx1p; ball1.Vy=Vy1p; ball1.Vz=Vz1p;
  ball2.Vx=Vx2p; ball2.Vy=Vy2p; ball2.Vz=Vz2p;
  //Integrar...
  //ball1.z+=ball1.Vz*dt; ball2.z+=ball2.Vz*dt;
}
void Collider::CollideBallsAlt(Ball & ball1, Ball & ball2, Crandom ran, double dt)
{
  double phi = 2*M_PI*ran.r(), P = ran.r();
  double theta = acos(sqrt(P));

  double Vx1 = ball1.Vx, Vy1 = ball1.Vy, Vz1 = ball1.Vz;
  double Vx2 = ball2.Vx, Vy2 = ball2.Vy, Vz2 = ball2.Vz;

  double Vcmx = (Vx1+Vx2)*0.5, Vcmy = (Vy1+Vy2)*0.5, Vcmz = (Vz1+Vz2)*0.5;
  double Ux1 = Vx1-Vcmx, Uy1 = Vy1-Vcmy, Uz1 = Vz1-Vcmz;
  double U1 = sqrt(Ux1*Ux1+Uy1*Uy1+Uz1*Uz1);
  double Ux2 = Vx2-Vcmx, Uy2 = Vy2-Vcmy, Uz2 = Vz2-Vcmz;
  double U2 = sqrt(Ux2*Ux2+Uy2*Uy2+Uz2*Uz2);
    
  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);

  double Vx1p = -U1*nx+Vcmx, Vx2p = U2*nx+Vcmx;
  double Vy1p = -U1*ny+Vcmy, Vy2p = U2*ny+Vcmy;
  double Vz1p = -U1*nz+Vcmz, Vz2p = U2*nz+Vcmz;

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
  double zi = 0, ball_zi;
  //dz = balls[N-1].z/M;
  //int sum_n = 0;   
  for(int l = 0;l<M;l++)
    {
      zi = l*dz;
      for(int i = 0; i<N; i++)
	{
	  ball_zi = balls[i].z;
	  if(ball_zi >= zi && ball_zi < zi+dz)
	    n[l]++;
	}
      
      /*sum_n += n[l];
      if(sum_n >= N)
      break;*/
    }
  
}
void Collider::VelocityDistribution(Ball * balls)
{
  int j = 0;
  double vi = 0, ball_vi=0;
  dv = 1e-3;
  //int sum_n = 0;   
  for(int l = 0;l<Mv;l++)
    {
      vi = (l-Mv/2+1)*dv;
      for(int i = 0; i<N; i++)
	{
	  ball_vi = balls[i].Vz;
	  if(ball_vi >= vi && ball_vi < vi+dv)
	    f[l]+=1./N;
	}
      //cout << vi << "\t" << f[l] << "\n";
    }
  
}
void Collider::PrintDensity(void)
{
  double zi = 0;
  for(int l = 0;l<M;l++)
    {
      zi = l*dz;
      cout << zi << "\t" << n[l] << "\n";
    }
}
void Collider::PrintVelocities(void)
{
  double vi = 0;
  for(int l = 0;l<Mv;l++)
    {
      vi = (l-Mv/2+1)*dz;
      cout << vi << "\t" << f[l] << "\n";
    }
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
double Collider::TotalEnergy(Ball * balls)
{
  double K=0.,U=0.,E=0.;
    for(int i=0;i<=N;i++)
      {
	K += 0.5*(balls[i].Vx*balls[i].Vx + balls[i].Vy*balls[i].Vy + balls[i].Vz*balls[i].Vz);
	U += g*balls[i].z;
      }
    E = K + U;
    return E;
}
int main()
{cout << fixed;cout.precision(8);
  Crandom rand64(1);
  Ball balls[N];
  Collider colls;

  double tColl = 0, time = 0; 
  int indexColl = 0;
  int collisions = 0;

  double E0 = 0.;
  bool areBallsSorted = true;

  //1. Inicializar las posiciones de las particulas.
  colls.InitPositions(balls, rand64);
  //InicieAnimacion();
  //cout << "coll" << "\t" << "time" << "\t" << "z1" << "\t" <<  "z2" << "\t" <<  "z3" << "\t" <<  "z4" << "\t" <<  "z5" << "\n" ;
  while(time <= 10000*tf)
    {
      // 0. Imprimir.
      //cout << collisions << "\t" << time;
      /*for(int i=0;i<N;i+=2)
	{
	  cout << "\t" <<  balls[i].GetZ();
	}
	cout << "\n" ;*/

      //2. Calcular todos los tiempos de colision entre i e i-1.
      colls.CollisionTimeBound(balls[0]);
      for(int i = 1;i<N;i++)
	{
	  colls.GetCollisionTime(balls[i-1],balls[i]);
	}
      //3. Obtener el tmin e I de colision.
      tColl = colls.GetTmin();
      if(tColl == 0)
	{
	  cout << "collision time is zero!\n";
	  break;
	}
      if(tColl >= 1e10)
	{
	  cout << "collision time is infinite!\n";
	  break;
	}
      indexColl = colls.GetIndex();
      //4. Integrar las particulas con dt=tmin.
      for(int i = 0;i<N;i++)
	{
	  balls[i].MovementEq();
	  balls[i].Integrate(tColl);
	}
      //5. Realizar la colision entre I e I-1.
      if(indexColl == 0){
	colls.CollideBound(balls[indexColl]);
	//balls[indexColl].Integrate(tColl);
      }
      else{
	colls.CollideBalls(balls[indexColl-1], balls[indexColl], rand64, tColl);
	//balls[indexColl-1].Integrate(tColl);
	//balls[indexColl].Integrate(tColl);
      }
      //6. Imprimir y actualizar el siguiente paso.
      /*InicieCuadro();
      for(int i=0;i<N;i++)
	  balls[i].PrintBall();
	  TermineCuadro();*/
      if(collisions > 5000 && collisions % 50 == 0)
	colls.MeasureDensity(balls);
      time += tColl;
      collisions++;
      //E0 = colls.TotalEnergy(balls);

    }
  //colls.MeasureDensity(balls);
  colls.PrintDensity();
  //for(int i=0;i<N;i++)cout << balls[i].GetI() << "\t" << balls[i].GetZ() << "\n";
  //cout << collisions << endl;
  return 0;
}
