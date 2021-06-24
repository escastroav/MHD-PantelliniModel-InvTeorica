#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int N = 50;
const double R0 = 1;

const double g = 0.98*R0;
const double H = N*R0*0.2; // kT/mg=rmax/5

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
  cout
    << z << "\t"
    << Vx << "\t" 
    << Vy << "\t" 
    << Vz << "\n"; 
}
class Collider
{
private:
  double t[N];
  double tmin=0.0;
  double I = 0;
public:
  Collider(void);
  double GetTmin(void){return tmin;};
  int GetIndex(void){return I;};
  void GetCollisionTime(Ball & ball1, Ball & ball2);
  void CollisionTimeBound(Ball & ball0);
  void CollideBalls(Ball & ball1, Ball & ball2, Crandom rand);
  void CollideBound(Ball & ball0);
  void InitPositions(Ball * balls, Crandom ran);
};
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
      if(t12 < tmin)
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
void Collider::CollideBalls(Ball & ball1, Ball & ball2, Crandom ran)
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
}
void Collider::CollideBound(Ball & ball0)
{
  double Vz0 = (-1)*ball0.Vz;
  ball0.Vz=Vz0;
}
void Collider::InitPositions(Ball * balls, Crandom ran)
{
  double zi=R0*ran.r();
  double vix = 0, viy = 0, viz = 0;
  for(int i = 0; i<N; i++)
    {
      zi += R0*ran.r();
      //viz = ran.r();
      balls[i].InitBall(i, zi, vix, viy, viz);
    }
}
int main()
{cout.precision(8);
  Crandom rand64(1);
  Ball balls[N];
  Collider colls;

  double tColl = 0, time = 0; 
  int indexColl = 0;

  //1. Inicializar las posiciones de las particulas.
  colls.InitPositions(balls, rand64);
  
  while(time < tf)
    {
      cout << time << "\t" << indexColl << "\t" << tColl << "\t" << balls[indexColl].GetZ() << "\t" << balls[indexColl].GetVz() << "\n" ;
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
	  cout << "tiempo de colision es cero!\n";
	  break;
	}
      if(tColl >= 1e10)
	{
	  cout << "tiempo de colision infinito!\n";
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
      if(indexColl == 0)
	colls.CollideBound(balls[indexColl]);
      else
	colls.CollideBalls(balls[indexColl], balls[indexColl-1], rand64);
      //6. Imprimir y actualizar el siguiente paso.
      /*
      for(int i=0;i<N;i++)
	{
	  balls[i].PrintBall();
	  };
	  time += tColl;*/
    }
  cout << tColl << endl;
  return 0;
}
