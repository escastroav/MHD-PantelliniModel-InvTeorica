#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#define _USE_MATH_DEFINES

using namespace std;
//Definición de constantes

const double size = 1;
const int N = 5;//*size;
const double R0 = 1/size;
const int seed = 0;
const double h = N*R0;
const double g = 3*R0;
const double H = h*0.2; // kT/mg=rmax/5
const double tf = sqrt(H/g);
const double eps = 1.0*pow(10, -10); //Para usar cuando comparamos doubles en los for

//generación de una variable que tenga los atributos de una partícula
struct particle{ 
  int id = 0;
  double z = 0.0;
  double Vx = 0.0;
  double Vy = 0.0;
  double Vz = 0.0;
  double az = -g;
};

struct collision{
  double tmin = 0.0;
  int I = 0; //Identificará la partícula que colisiona en el tiempo tmin con ii+1
};
  
void initial_values_particles (vector<particle> & particles);
void movement_equations (vector<particle> & particles, double dt);
void sort_particles (vector<particle> & particles);
void Get_Collision_Time(vector<particle> & particles, vector<double> & times, collision & collisions);

int main ()
{
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);
  vector<particle> particles (N);
  initial_values_particles (particles);
  collision time_collisions;
  
  //Imprime primero las partículas con su información de datos dados aleatoriamente. Luego, evoluciona el sistema un tiempo dt y calcula los nuevos parámetros. Posteriormente, mira la información de las alturas de cada partícula y las ordena, cambiandoles el id. Luego se añade los tiempos de colision y se escoge el mínimo (No se ha tenido en cuenta la colisión con las fronteras)
  
  for (int ii=0; ii<N; ii++)
    {
      std::cout<<ii<<"\t"<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }

  movement_equations (particles, 1.2);

  for (int ii=0; ii<N; ii++)
    {
      std::cout<<ii<<"\t"<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }
  
  sort_particles (particles);
  
  for(int ii=0; ii<N; ii++)
    {
      std::cout<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }

  vector<double> times (N-1, 0.0);

  Get_Collision_Time(particles, times, time_collisions);
  
  for (int ii=0; ii<(N-1); ii++) //Para evitar un segmentation fault
    {
      std::cout<<ii<<"\t"<<times[ii]<<std::endl;
    }

  std::cout<<time_collisions.I<<"\t"<<time_collisions.tmin<<"\t"<<*min_element(times.begin(), times.end())<<std::endl;
  
  return 0;
}

void initial_values_particles (vector<particle> & particles) //Tiene como objetivo dar valores iniciales aleatorios a las partículas
{
  //Enumera las partículas
  for (int ii=0; ii<N; ii++)
    {
      particles[ii].id = ii;
      //genera numeros distribuidos uniformemente (Toca editar porque los perfiles de velocidad tienen una distribución maxwelliana)
      std::mt19937 gen(ii);
      std::uniform_real_distribution<double> dis(0, 5.0);
      std::uniform_real_distribution<double> disp(0, 10.0); 
      particles[ii].z = disp(gen);
      particles[ii].Vx = dis(gen);
      particles[ii].Vy = dis(gen);
      particles[ii].Vz = dis(gen);
    }
}

void movement_equations (vector<particle> & particles, double dt)
{
  for (int ii=0; ii<N; ii++)
    {
      particles[ii].z += particles[ii].Vz*dt-0.5*g*pow(dt, 2);
      particles[ii].Vz += -g*dt;
    }
}

void sort_particles (vector<particle> & particles)
{
  //ordena de menor a mayor los elementos del vector 
  sort(particles.begin(), particles.end(), [](particle a, particle b) {
		return a.z < b.z;
	});
  
  //cambia el id de la partícula ordenada
  for (int ii=0; ii<N; ii++)
    {
      particles[ii].id = ii;
    }
}

void Get_Collision_Time(vector<particle> & particles, vector<double> & times, collision & collisions)
{
  //vector<double> times (N-1, 0.0); //Son N-1 posibles tiempos a calcular ya que van de manera consecutiva (3 partículas equivale al tiempo t12 y t23)
  
  for (int ii=0; ii<(N-1); ii++) //Pues el N-1 asegura que no haya un segmentation fault en el for
    {
      double vrel = (particles[ii+1].Vz-particles[ii].Vz);
      if(abs(vrel) > eps)
	{
	  times[ii] = abs((particles[ii+1].z-particles[ii].z)/vrel);
	}
      else{
        times[ii] = 1.0*pow(10, 10);
      }
    }
  // Hacer una estructura collide con info de la partícula I y el tiempo minimo tmin
  
  collisions.tmin = *min_element(times.begin(), times.end());
  
  for (int ii=0; ii<(N-1); ii++) //Podría verse si se puede optimizar ya que se hace un for solamente para calcular I pudiendose obtener en el for previo.
    {
      if (times[ii] == collisions.tmin)
	{
	  collisions.I = ii+1;  //Pues los tiempos son tomados como ii a I-1 y ii+1 como I ya que ii inicia en cero
	  break;
	}
    }
}

//Revisar notación
/*
void Collide_particles (vector<particle> particles)
{
  double phi = 2*M_PI*ran.r(), P = ran.r();//Número aleatorio entre 0 a 2pi
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
*/
