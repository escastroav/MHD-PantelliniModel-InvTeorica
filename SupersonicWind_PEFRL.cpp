#include <iostream>
#include <cmath>
#include "Random64.h"
#include <iomanip>
#include <cstdlib>
using namespace std;

const double size = 1;
const int N = 50*size;
const double R0 = 1./size;

const double h = N*R0;
const double g = 9.8e-3*R0;
const double H = h*.2; // kT/mg=rmax/5

const double tf = sqrt(H/g);

const double mp = 1.;
const double me = 1e-2;

const double e = 1e-3;

const double U_ab = 0.5;

const double Zi = 0.1786178958448091;
const double Lambda = -0.2123418310626054;
const double Xi = -0.06626458266981849;

const double Coeff1 = (1-2*Lambda)*0.5;
const double Coeff2 = (1-2*(Xi+Zi));

class Ball
{
private:
  int i = 0;
  double z=0,
    Vx=0,
    Vy=0,
    Vz=0;
  double az=0;
  double q=0;
  double m=0;
  double E=0;
public:
  void InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0,double q0, double m0, double E0);
  double GetZ(void){return z;};void SetZ(double zi){z=zi;};
  double GetQ(void){return q;};void SetQ(double qi){q=qi;};  
  double GetE(void){return E;};void SetE(double Ei){E=Ei;};  
  double GetVx(void){return Vx;};void SetVx(double Vxi){Vx=Vxi;};
  double GetVy(void){return Vy;};void SetVy(double Vyi){Vy=Vyi;};
  double GetVz(void){return Vz;};void SetVz(double Vzi){Vz=Vzi;};
  double GetAz(void){return az;};  
  int GetI(void){return i;};
  void MovementEq(void);
  void IntegrateZ(double dt, double coeff);
  void IntegrateV(double dt, double coeff);
  void PrintBall(void);
  friend class Collider;
};
void Ball::InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0,double q0, double m0, double E0)
{i=i0; z=z0; Vx=Vx0; Vy=Vy0; Vz=Vz0;q=q0;m=mp;E=E0;}
void Ball::MovementEq(void)
{
  az = -g+q*E/m;
}
void Ball::IntegrateZ(double dt, double coeff)
{
  z += Vz*dt*coeff;
}
void Ball::IntegrateV(double dt, double coeff)
{
  Vz += az*dt*coeff;
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
  double * ne = nullptr;
  double * np = nullptr;
  double dz = 0.245*H;
  int M = 0;
public:
  Collider(void);
  ~Collider(void);
  double GetTmin(void){return tmin;};
  int GetIndex(void){return I;};
  int GetLowIndex(Ball * balls);
  void ElectricField(Ball * balls, double epsilon, double E0);
  void GetCollisionTime(Ball & ball1, Ball & ball2);
  void CollisionTimeBound(Ball & ball0);
  void CollideBalls(Ball & ball1, Ball & ball2, Crandom rand, double dt);
  void CollideBound(Ball & ball0);
  void InitPositions(Ball * balls, Crandom ran, double E0);
  double ChargeNeutrality(Ball * balls);
  void MeasureDensity(Ball * balls);
  void PrintDensity(void);
  double TotalEnergy(Ball * balls);
  bool AreBallsSorted(Ball * balls);
  double CollisionProbability(Ball & ball1, Ball & ball2);
};
Collider::Collider(void)
{
  M = 100;
  dz = h / M;
  ne = new double[M];
  np = new double[M];
  for(int i=0;i<N;i++)
    t[i] = 0;
}
Collider::~Collider(void)
{
  delete ne; delete np;
}
void Collider::ElectricField(Ball * balls,double epsilon, double E0)
{
  balls[0].E = E0 + (balls[N-1].q-balls[0].q)*epsilon*0.5;
  for(int i=0;i<N;i++)
      balls[i+1].E = balls[i].E + balls[i].q*epsilon*0.5; 
    
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
  double theta = acos(sqrt(P));
  double m1 = ball1.m, m2 = ball2.m;

  double Vx1 = ball1.Vx, Vy1 = ball1.Vy, Vz1 = ball1.Vz;
  double Vx2 = ball2.Vx, Vy2 = ball2.Vy, Vz2 = ball2.Vz;

  double Ux = Vx2-Vx1, Uy = Vy2-Vy1, Uz = Vz2-Vz1;
  double  U = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);

  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);

  double Vx1p = (m1*Vx1+m2*Vx2+U*m2*nx)/(m1+m2), Vx2p = (m1*Vx1+m2*Vx2-U*m1*nx)/(m1+m2);
  double Vy1p = (m1*Vy1+m2*Vy2+U*m2*ny)/(m1+m2), Vy2p = (m1*Vy1+m2*Vy2-U*m1*ny)/(m1+m2);
  double Vz1p = (m1*Vz1+m2*Vz2+U*m2*nz)/(m1+m2), Vz2p = (m1*Vz1+m2*Vz2-U*m1*nz)/(m1+m2);

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
void Collider::InitPositions(Ball * balls, Crandom ran, double E0)
{
  double zi=R0, qi = 0, mi = 0;
  double vix = 0, viy = 0, viz = 0;
  for(int i = 0; i<N; i++)
    {
      zi += R0;
      if(i%2==0){
	qi = e; mi = mp;}
      else{
	qi = -e; mi = me;}
      balls[i].InitBall(i, zi, vix, viy, viz, qi, mi, E0);
    }
}
int Collider::GetLowIndex(Ball * balls)
{
  double zmin = balls[0].z;
  int imin = balls[0].i;
  for(int i=0;i<N;i++)
    {
      if(balls[i].z < zmin)
	{
	  zmin = balls[i].z;
	  imin = balls[i].i;
	}
    }
  return imin;
}
double Collider::ChargeNeutrality(Ball * balls)
{
  double moments = 0;
  for(int i=0;i<N;i++)
    {
      moments += balls[i].q*balls[i].z; 
    }
  return moments/h;
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
	    {
	      if(balls[i].q < 0)
		ne[l]+=1;
	      else
		np[l]+=1;
	    }
	}
      
      /*sum_n += n[l];
      if(sum_n >= N)
      break;*/
    }
  
}
void Collider::PrintDensity(void)
{
  double zi = 0;
  for(int l = 0;l<M;l++)
    {
      zi = l*dz;
      cout << zi << "\t" << ne[l] << "\t" << np[l] << "\n";
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
double Collider::CollisionProbability(Ball & ball1, Ball & ball2)
{
  double ux = ball2.Vx - ball1.Vx, uy = ball2.Vy - ball1.Vy, uz = ball2.Vz - ball1.Vz;
  double u = sqrt(ux*ux+uy*uy+uz*uz);
  double R_ab = U_ab/u;  
  if(u < U_ab)
    return 1.;
  else
    return (R_ab*R_ab*R_ab*R_ab);
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
  Crandom collRand(5);
  Ball balls[N];
  Collider colls;

  double tColl = 0, time = 0; 
  int indexColl = 0;
  int collisions = 0;

  int imin=0;

  double E0 = 1., eps = 1., polar = 0;
  double Pu = 0, chanceColl = 0;
  bool areBallsSorted = true;

  //1. Inicializar las posiciones de las particulas.
  colls.InitPositions(balls, rand64, E0);
  //InicieAnimacion();
  //cout << "coll" << "\t" << "time" <<  "\t" << "t_min" << "\t" << "z1" << "\t" <<  "z2" << "\t" <<  "z3" << "\t" <<  "z4" << "\t" <<  "z5" << "\n" ;
  //double tInitial = sqrt(2*balls[0].GetZ()/g);
  //for(int i = 0;i<N;i++) balls[i].StartPosition(tInitial);
  while(collisions <= 100000)
    {
      //0. Imprimir.
      //E0 = colls.TotalEnergy(balls);
      /*
	for(int i=0;i<N;i+=2)
	  {
	    cout << time << "\t" << balls[i].GetZ();
	  }
	
	cout  << "\n";
      */
	//<< "\t" << balls[2].GetZ()
	//<< "\t" << balls[3].GetZ()
	/*	   << "\t" << balls[4].GetZ()
	   << "\t" << balls[5].GetZ()
	   << "\t" << balls[6].GetZ()
	   << "\t" << balls[7].GetZ()
	   << "\t" << balls[8].GetZ()
	   << "\t" << balls[9].GetZ() <<"\n" ;*/
      //2. Calcular todos los tiempos de colision entre i e i-1.
      imin = colls.GetLowIndex(balls);
      colls.CollisionTimeBound(balls[imin]);
      for(int i = 0;i<N;i++)
	{
	  colls.GetCollisionTime(balls[i-1],balls[i]);
	}
      //3. Obtener el tmin e I de colision.
      tColl = colls.GetTmin();indexColl = colls.GetIndex();
      if(tColl == 0)
	{
	  cout << "collision time is zero!\t" << indexColl << "\t" << balls[indexColl].GetZ() << "\n";
	  break;
	}
      if(tColl >= 1e10)
	{
	  cout << "collision time is infinite!\t" << indexColl << "\t" << balls[indexColl].GetZ() << "\n" ;
	  break;
	}
      
      
      //4. Integrar las particulas con dt=tmin.
      colls.ElectricField(balls,eps,E0);
      for(int i = 0;i<N;i++)
	{
	  balls[i].IntegrateZ(tColl,Zi);
	  balls[i].MovementEq();
	  balls[i].IntegrateV(tColl,Coeff1);
	  balls[i].IntegrateZ(tColl,Xi);
	  balls[i].MovementEq();
	  balls[i].IntegrateV(tColl,Lambda);
	  balls[i].IntegrateZ(tColl,Coeff2);
	  balls[i].MovementEq();
	  balls[i].IntegrateV(tColl,Lambda);
	  balls[i].IntegrateZ(tColl,Xi);
	  balls[i].MovementEq();
	  balls[i].IntegrateV(tColl,Coeff1);
	  balls[i].IntegrateZ(tColl,Zi);
	}
      polar = colls.ChargeNeutrality(balls);
      //5. Realizar la colision entre I e I-1.
      if(indexColl == 0){
	colls.CollideBound(balls[indexColl]);
	//balls[indexColl].Integrate(tColl);
      }
      else{
	//srand(indexColl);//collRand.Reset((unsigned long long)rand());
        //chanceColl = rand()/RAND_MAX;
	//Pu = 1;//colls.CollisionProbability(balls[indexColl], balls[indexColl-1]);
	//if(chanceColl < Pu)
	  colls.CollideBalls(balls[indexColl], balls[indexColl-1], rand64, tColl);
	//balls[indexColl-1].Integrate(tColl);
	//balls[indexColl].Integrate(tColl);
      }
      //6. Imprimir y actualizar el siguiente paso.
      /*InicieCuadro();
      for(int i=0;i<N;i++)
	  balls[i].PrintBall();
	  TermineCuadro();*/
      /*if(collisions > 5000 && collisions % 50 == 0)
	colls.MeasureDensity(balls);*/
	time += tColl;
      collisions++;
      
    }
  //colls.PrintDensity();
  for(int i=0;i<N;i++)cout << balls[i].GetI() << "\t" << balls[i].GetQ() << "\t" << balls[i].GetZ() << "\t" << balls[i].GetAz() << "\n";
  cout << collisions << endl;
  return 0;
}
