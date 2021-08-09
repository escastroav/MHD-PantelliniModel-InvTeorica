#include <iostream>
#include <cmath>
#include "Random64.h"
#include <iomanip>
#include <cstdlib>
using namespace std;

const double gmma = 23;

const double Rsun = 1.;//1 unidad de longitud (u.a.L)= Rsun(SI)/10 = 6.9634e7 m
const int N = 100;//80 protones y 80 electrones. N no sebe ser impar

const double L = 5*Rsun;
const double R0 = L/N;

const double g = 3.934e-6;// u.a.L/s^2

const double tf = sqrt(L/g);//504.12 s

const double mp = 1.6e-8;//no tiene que ser el real pero ajá 
const double me = 1e-2*mp;//mass ratio me/mp = 100

const double q0 = 1.6e-3;//??

const double T0 = 5e5;// Kelvin
const double TL = 6.4e5;//

const double kB = g*Rsun*mp/(2*gmma*T0);//Ec 12 del paper On temperature Profile...

const double vm = sqrt(2*kB*T0/me);
const double U_ab = 1e-1;

const double ET = g*mp*(1-1e-2)/(2*q0);

class Ball
{
private:
  int i = 0;
  double z=0,
    Vx=0,
    Vy=0,
    Vz=0;
  double az=0;
  double m = mp;
  double q = q0;

  double E = 0;
public:
  void InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0,double qs, double ms);
  double GetZ(void){return z;};void SetZ(double zi){z=zi;};
  double GetM(void){return m;};void SetM(double mi){m=mi;};
  double GetQ(void){return q;};void SetQ(double qi){q=qi;};
  double GetE(void){return E;};void SetE(double Ei){E=Ei;};
  double GetVx(void){return Vx;};void SetVx(double Vxi){Vx=Vxi;};
  double GetVy(void){return Vy;};void SetVy(double Vyi){Vy=Vyi;};
  double GetVz(void){return Vz;};void SetVz(double Vzi){Vz=Vzi;};
  int GetI(void){return i;};
  void MovementEq(double epsilon, double E0);
  void ElectricField(double epsilon, double E0);
  void Integrate(double dt);
  friend class Collider;
};
void Ball::InitBall(int i0, double z0, double Vx0, double Vy0, double Vz0, double qs, double ms)
{
  i=i0; z=z0; Vx=Vx0; Vy=Vy0; Vz=Vz0; q=qs; m=ms;
  E = ET;
  az = -g + q*E/m;
}
void Ball::ElectricField(double epsilon, double E0)
{
  E += q*epsilon*0.5 + E0;
}
void Ball::MovementEq(double epsilon,double E0)
{
  ElectricField(epsilon,E0);
  az = -g + q*E/m;
}
void Ball::Integrate(double dt)
{
  z += Vz*dt + az*dt*dt*0.5;
  Vz += az*dt;
}
class Collider
{
private:
  double t[N+1];
  double E[N];
  double tmin=0.0;
  int I = 0;
  double * n = nullptr;
  double dz = 0.245*L;
  int M = 0;
public:
  Collider(void);
  ~Collider(void);
  double GetTmin(void){return tmin;};
  int GetIndex(void){return I;};
  int LowestParticleIndex(Ball * balls);
  int HighestParticleIndex(Ball * balls);
  void GetCollisionTime(Ball & ball1, Ball & ball2);
  void GetElectricField(Ball * ball, double epsilon);
  void CollisionTimeGround(Ball & ball0);
  void CollisionTimeCeil(Ball & ball0);
  void CollideBalls(Ball & ball1, Ball & ball2, Crandom rand, double dt);
  void CollideGround(Ball & ball0, Crandom ran);
  void CollideCeil(Ball & ball0, Crandom ran);
  void InitPositions(Ball * balls, Crandom ran);
  void MeasureDensity(Ball * balls);
  void PrintDensity(void);
  void SortParticles(Ball * balls);
  void ShowTimes(void);
  double CollisionProbability(Ball & ball1, Ball & ball2);
  double TotalEnergy(Ball * balls);
  double ChargeNeutrality(Ball * balls);

};
Collider::Collider(void)
{
  M = 100;
  dz = L / M;
  n = new double[M];
  for(int i=0;i<N;i++)
    t[i] = 0;
}
Collider::~Collider(void)
{
  delete n;
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
void Collider::GetElectricField(Ball* balls,double epsilon)
{
  double Q0 = balls[0].q, QN = balls[N-1].q;
  E[0] = ET + (QN-Q0)*epsilon*0.5;
  for(int i=1;i<N;i++)
    {
      E[i] = E[i-1] + balls[i].q*epsilon;
      balls[i].E = E[i];
    }
}
void Collider::CollisionTimeGround(Ball & ball0)
{
  double t0=0, z0 = ball0.z, Vz0 = ball0.Vz, az0 = ball0.az;
  double disc = Vz0*Vz0-2*az0*z0, t1=0,t2=0;
  /*if(abs(z0) < 1e-10)
    t0 = 1e10;*/
  if(az0 == 0)
    {
      if(Vz0<0)
	t0 = -z0/Vz0;
      else
	t0 = 1e10;
    }     
  else if(disc >= 0)
    {
      t1 = (- Vz0 - sqrt(disc))/az0;
      t2 = (- Vz0 + sqrt(disc))/az0;
      if(t1*t2 < 0)
	t0 = (t1>0)? t1 : t2;
      else if(t1 > 0)
	t0 = (t1<t2)? t1 : t2;
      else
	t0 = 1e10;
    }
  else
    t0 = 1e10;

  if(t0>0 && t0 != 1e10)
    {tmin = t0;
      I = ball0.i;}
  t[0] = t0;
}
void Collider::CollisionTimeCeil(Ball & ball0)
{
  double t0=0, z0 = ball0.z-L, Vz0 = ball0.Vz, az0 = ball0.az;
  double disc = Vz0*Vz0-2*az0*z0, t1=0,t2=0;
  /*if(abs(z0) < 1e-10)
    t0 = 1e10;*/
  if(az0 == 0)
    {
      if(Vz0>0)
	t0 = -z0/Vz0;
      else
	t0 = 1e10;
    }     
  else if(disc >= 0)
    {
      t1 = (- Vz0 - sqrt(disc))/az0;
      t2 = (- Vz0 + sqrt(disc))/az0;
      if(t1*t2 < 0)
	t0 = (t1>0)? t1 : t2;
      else if(t1 > 0)
	t0 = (t1<t2)? t1 : t2;
      else
	t0 = 1e10;
    }
  else
    t0 = 1e10;
  if(t0>0 && t0<tmin)    
    {tmin = t0;  
      I = ball0.i+1;}
  t[N] = t0;
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

  double Vx1p = (m1*Vx1+m2*Vx2+m2*U*nx)/(m1+m2), Vx2p = (m1*Vx1+m2*Vx2-m1*U*nx)/(m1+m2);
  double Vy1p = (m1*Vy1+m2*Vy2+m2*U*ny)/(m1+m2), Vy2p = (m1*Vy1+m2*Vy2-m1*U*ny)/(m1+m2);
  double Vz1p = (m1*Vz1+m2*Vz2+m2*U*nz)/(m1+m2), Vz2p = (m1*Vz1+m2*Vz2-m1*U*nz)/(m1+m2);

  ball1.Vx=Vx1p; ball1.Vy=Vy1p; ball1.Vz=Vz1p;
  ball2.Vx=Vx2p; ball2.Vy=Vy2p; ball2.Vz=Vz2p;
  //Integrar...
  //ball1.z+=ball1.Vz*dt; ball2.z+=ball2.Vz*dt;
}
void Collider::CollideGround(Ball & ball0, Crandom ran)//Maxwellian distribution
{
  double alpha = ball0.m/(2*kB*T0), sigma = sqrt(1./(2*alpha));
  ran.Reset((unsigned long long)rand());//Actualziar seed para que ran.r() no sea el mismo.
  double V0 = abs(ran.gauss(0.,sigma));
  ran.Reset((unsigned long long)rand());
  double P_th=ran.r(),theta=acos(sqrt(P_th));
  /* No es necesario este ciclo... 
while(true)//Intento de generar variable aleatorio con distribución sin^2(theta)
    {//Alguien que se acuerde de sus clases de probabilidad, pa ver como se hace esto!!!
      ran.Reset((unsigned long long)rand());
      P_th = asin(sqrt(ran.r()));
      ran.Reset((unsigned long long)rand());
      test_th = ran.r();
      if(test_th < P_th)
	{
	  theta = test_th;
	  break;
	}
	}*/
  
  ran.Reset((unsigned long long)rand());
  double phi = 2*M_PI*ran.r();
  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);
  ball0.z = 0;
  ball0.Vx=V0*nx;
  ball0.Vy=V0*ny;
  ball0.Vz=V0*nz;
}
void Collider::CollideCeil(Ball & ball0, Crandom ran)//Maxwellian distribution
{
  double alpha = ball0.m/(2*kB*T0), sigma = sqrt(1./(2*alpha));
  ran.Reset((unsigned long long)rand());//Actualziar seed para que ran.r() no sea el mismo.
  double V0 = abs(ran.gauss(0.,sigma));
  ran.Reset((unsigned long long)rand());
  double P_th=ran.r(),theta=acos(sqrt(P_th));
  /* No es necesario este ciclo... 
while(true)//Intento de generar variable aleatorio con distribución sin^2(theta)
    {//Alguien que se acuerde de sus clases de probabilidad, pa ver como se hace esto!!!
      ran.Reset((unsigned long long)rand());
      P_th = asin(sqrt(ran.r()));
      ran.Reset((unsigned long long)rand());
      test_th = ran.r();
      if(test_th < P_th)
	{
	  theta = test_th;
	  break;
	}
	}*/
  
  ran.Reset((unsigned long long)rand());
  double phi = 2*M_PI*ran.r();
  double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);
  ball0.z=L; 
  ball0.Vx=V0*nx;
  ball0.Vy=V0*ny;
  ball0.Vz=-V0*nz;
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
void Collider::InitPositions(Ball * balls, Crandom ran)
{
  double zi=R0;
  double chance=0;
  double qs=0, ms=0;
  double vix=0,viy=0,viz=0;
  for(int i = 0; i<N; i++)
    {
      if(i % 2 == 0){//es proton
	qs = q0;
	ms = mp;
      }else{//es electron
	qs = -q0;
	ms = me;
      }
      balls[i].InitBall(i, zi, vix, viy, viz, qs, ms);
      ran.Reset((unsigned long long)rand());
      zi += R0*ran.r();
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
	    n[l]+=1;
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
      cout << zi << "\t" << n[l] << "\n";
    }
}
double Collider::ChargeNeutrality(Ball * balls)
{
  double moments = 0;
  for(int i=0;i<N;i++)
    {
      moments += balls[i].q*balls[i].z; 
    }
  return moments/L;
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
void Collider::SortParticles(Ball * balls)
{
  Ball ballSwap = balls[0];
  for(int j=0;j<N-1;j++)
    {
      if(balls[j+1].z < balls[j].z)
	{
	  ballSwap = balls[j+1];
	  balls[j+1] = balls[j];
	  balls[j] = ballSwap;
	}
    }
}
int Collider::LowestParticleIndex(Ball * balls)
{
  double zmin = balls[0].z;
  double I_min = balls[0].i;
  for(int i=1;i<N;i++)
    {
      if(balls[i].z < zmin){
	zmin = balls[i].z;
	I_min = balls[i].i;
      }
    }
  return I_min;
}
int Collider::HighestParticleIndex(Ball * balls)
{
  double zmax = balls[N-1].z;
  double I_max = balls[N-1].i;
  for(int i=N-2;i>=0;i--)
    {
      if(balls[i].z > zmax){
	zmax = balls[i].z;
	I_max = balls[i].i;
      }
    }
  return I_max;
}
void Collider::ShowTimes(void)
{
  cout << "times" << endl;
  for(int i=0;i<N+1;i++)
    cout << t[i] << "\t";
  cout << endl;
}
int main()
{cout << fixed;cout.precision(8);
  Crandom rand64(1);
  Ball balls[N];
  Collider colls;

  int imin=0,imax=0;
  double tColl = 0, time = 0; 
  int indexColl = 0;
  int collisions = 0;

  double E0 = ET;
  double eps = mp*vm*vm/(R0*q0*q0*10);//...

  double collP = 0, collR = 0;
  //1. Inicializar las posiciones de las particulas.
  colls.InitPositions(balls, rand64);
	
  double polarization=0,minPol=q0*L;//??
  /*cout << "index" <<"\t\t"<< "m" <<"\t\t"<< "q" <<"\t\t"<< "z" <<"\t\t"<< "vz" << "\t\t" << "Ei" << "\n"; 
  for(int i=0;i<N;i++)
    cout << balls[i].GetI() << "\t"
	 << balls[i].GetM() << "\t"
	 << balls[i].GetQ() << "\t"
	 << balls[i].GetZ() << "\t"
	 << balls[i].GetVz() << "\t"
	 << balls[i].GetE() << "\n"; 
  cout << "-------------------------------------------------" << endl;
  */
  while(collisions <= 15000)
    {
      imin = colls.LowestParticleIndex(balls);
      colls.CollisionTimeGround(balls[imin]);
      for(int i = 1;i<N;i++)
	{
	  colls.GetCollisionTime(balls[i-1],balls[i]);
	}
      imax = colls.HighestParticleIndex(balls);
      colls.CollisionTimeCeil(balls[imax]);
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
      if(tColl < 0)
	{
	  cout << "collision time is negative!\n";
	  break;
	}
      indexColl = colls.GetIndex();
      /*colls.ShowTimes();
      cout << "collison" << "\t\t" << "tColl" << "\t\t" << "indexColl" << endl;
      cout << collisions << "\t" << tColl << "\t" << indexColl << endl;*/
      //4. Integrar las particulas con dt=tmin.
      for(int i = 0;i<N;i++){
	  balls[i].Integrate(tColl);
	  balls[i].MovementEq(eps,E0);}
      //5. Realizar la colision entre I e I-1. 
      if(indexColl == 0){
	colls.CollideGround(balls[indexColl],rand64);
      }
      else if(indexColl == N){
	colls.CollideCeil(balls[indexColl-1],rand64);
      }
      else{
	rand64.Reset((unsigned long long)rand());
	collP = colls.CollisionProbability(balls[indexColl], balls[indexColl-1]);
	collR = rand64.r();
	if(collR < collP)//decidir si colisiona o no. (no estoy seguro si es asi pero tiene sentido.)
	  colls.CollideBalls(balls[indexColl], balls[indexColl-1], rand64, tColl);
	/*else
	  cout << "particles didn't collided!\n";*/
      }
      //6. Verificar si el sistema es neutro electricamente
      polarization = colls.ChargeNeutrality(balls);
      if(polarization > minPol){
	cout << "polarization is high! : " << polarization << "\n";
	E0 -= polarization;}
      if(polarization < -minPol){
	cout << "polarization is low! : " << polarization << "\n";
	E0 -= polarization;}
      
      colls.SortParticles(balls);
      
      //7. Imprimir y actualizar el siguiente paso.
      //cout << "index" <<"\t\t"<< "m" <<"\t\t"<< "q" <<"\t\t"<< "z" <<"\t\t"<< "vz" << "\t\t" << "Ei" <<"\n"; 
      /*for(int i=0;i<N;i++)
	cout << balls[i].GetI() << "\t"
	     << balls[i].GetM() << "\t"
	     << balls[i].GetQ() << "\t"
	     << balls[i].GetZ() << "\t"
	     << balls[i].GetE() << "\t"
	     << balls[i].GetVz() << "\n";
      cout << "-------------------------------------------------" << endl;
      //if(collisions > 5000 && collisions % 50 == 0)
      //	colls.MeasureDensity(balls);*/
      time += tColl;
      collisions++;
      
    }
  //colls.PrintDensity();
  cout << "index" <<"\t\t"<< "m" <<"\t\t"<< "q" <<"\t\t"<< "z" <<"\t\t"<< "vz" << "\t\t" << "Ei" <<"\n"; 
      for(int i=0;i<N;i++)
	cout << balls[i].GetI() << "\t"
	     << balls[i].GetM() << "\t"
	     << balls[i].GetQ() << "\t"
	     << balls[i].GetZ() << "\t"
	     << balls[i].GetE()-ET << "\t"
	     << balls[i].GetVz() << "\n";
      cout << "-------------------------------------------------" << endl;
  return 0;
}
