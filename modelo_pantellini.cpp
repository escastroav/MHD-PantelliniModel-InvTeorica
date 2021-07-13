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

//generación de una variable que tenga los atributos de una partícula
struct particle{ 
  int id = 0;
  double z = 0.0;
  double Vx = 0.0;
  double Vy = 0.0;
  double Vz = 0.0;
  double az = -g;
};

void initial_values_particles (vector<particle> & particles);
void movement_equations (vector<particle> & particles, double dt);
void sort_particles (vector<particle> & particles);


int main ()
{
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);
  vector<particle> particles (N);
  initial_values_particles (particles);
  
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
  
  /*
  vector<double> height (N, 0.0);
  sort_particles (particles, height);
  for(int ii=0; ii<N; ii++)
    {
      std::cout<<height[ii]<<"\t";
    }
  */
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
