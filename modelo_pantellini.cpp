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
const double R0 = 5.0;
const double Rf = 50*R0;
const int seed = 0;
const double h = N*R0;
const double g = 3*R0;
const double H = h*0.2; // kT/mg=rmax/5
const double tf = sqrt(H/g);
const double eps = 1.0*pow(10, -10); //Para usar cuando comparamos doubles en los for
const double infty = 1.0*pow(10, 10); //Para expresar que los tiempos son muy grandes

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
  double tmin = 0.0; //Tiempo mínimo de colision entre partículas continuas
  double t_intb = 0.0; //Tiempo de colisión de la primera partícula con la forntera base.
  double t_extb = 0.0; //Tiempo de colisión con la frontera exterior de la ultima partícula.
  int I = 0; //Identificará la partícula que colisiona en el tiempo tmin con ii+1
};
  
void initial_values_particles (vector<particle> & particles);
void movement_equations (vector<particle> & particles, double dt);
void sort_particles (vector<particle> & particles);
void Get_Collision_Time(vector<particle> & particles, vector<double> & times, collision & collisions);
void collision_time_boundary (vector<particle> & particles, vector<double> & times, collision & collisions);
double Collide_particles_and_evolve_system (vector<particle> & particles, vector<double> & times, collision & collisions);

int main ()
{
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);
  vector<particle> particles (N);
  initial_values_particles (particles);
  collision time_collisions;
  
  //Imprime primero las partículas con su información de datos dados aleatoriamente. Luego, evoluciona el sistema un tiempo dt y calcula los nuevos parámetros. Posteriormente, mira la información de las alturas de cada partícula y las ordena, cambiandoles el id. Luego se añade los tiempos de colision y se escoge el mínimo (No se ha tenido en cuenta la colisión con las fronteras)

  std::cout<<"Se imprimen para 5 partículas sus valores asigandos aleatoriamente para id, z, Vx, Vy, Vz"<<std::endl;
  
  for (int ii=0; ii<N; ii++)
    {
      std::cout<<ii<<"\t"<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }

  movement_equations (particles, 0.5);

  std::cout<<"Se evoluciona el sistema un dt=0.5 y se imprimen para las mismas 5 partículas sus valores para id, z, Vx, Vy, Vz"<<std::endl;
  
  for (int ii=0; ii<N; ii++)
    {
      std::cout<<ii<<"\t"<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }
  
  sort_particles (particles);

  std::cout<<"Se organizan las partículas de acuerdo a su posición y se imprimen id, z, Vx, Vy, Vz"<<std::endl;
  
  for(int ii=0; ii<N; ii++)
    {
      std::cout<<particles[ii].id<<"\t"<<particles[ii].z<<"\t"<<particles[ii].Vx<<"\t"<<particles[ii].Vy<<"\t"<<particles[ii].Vz<<std::endl;
    }

  vector<double> times (N+1, 0.0); //En el espacio [N] y [N-1] (último dos espacios) guardará los tiempos de colisión con las fronteras int y ext.

  Get_Collision_Time(particles, times, time_collisions);
  collision_time_boundary (particles, times, time_collisions);

  std::cout<<"Se imprimen los tiempos de colisión (los primeros N-1 datos son colision entre particulas y los dos últimos son con fronteras)"<<std::endl;
  
  for (int ii=0; ii<(N+1); ii++) 
    {
      std::cout<<ii<<"\t"<<times[ii]<<std::endl;
    }

  std::cout<<"Se imprimen los tiempos mínimos de colisión (entre partículas y fronteras)"<<std::endl;
  
  std::cout<<time_collisions.I<<"\t"<<time_collisions.tmin<<"\t"<<time_collisions.t_intb<<"\t"<<time_collisions.t_extb<<"\t"<<*min_element(times.begin(), times.end())<<std::endl;
  
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
  
  for (int ii=0; ii<(N-1); ii++) //Pues el N-1 asegura que no se guarde ningun tiempo en el último espacio (designado para tiempos con las fronteras)
    {
      double vrel = (particles[ii+1].Vz-particles[ii].Vz);
      if(abs(vrel) > eps)
	{
	  times[ii] = abs((particles[ii+1].z-particles[ii].z)/vrel);
	}
      else{
        times[ii] = infty;
      }
    }
  // Hacer una estructura collide con info de la partícula I y el tiempo minimo tmin
  
  collisions.tmin = *min_element(times.begin(), times.end()-2); //No tiene en cuenta los dos últimos elementos porque esos tiempos están asignados para las fronteras
  
  for (int ii=0; ii<(N-1); ii++) //Podría verse si se puede optimizar ya que se hace un for solamente para calcular I pudiendose obtener en el for previo.
    {
      if (abs(times[ii]-collisions.tmin) < eps) //pues se comparan double entonces usar == podría fallar
	{
	  collisions.I = ii+1;  //Pues los tiempos son tomados como ii a I-1 y ii+1 como I ya que ii inicia en cero
	  break;
	}
    }
}

void collision_time_boundary (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  //Frontera interna

  double t_base = 0.0, Vz0 = particles[0].Vz, z0 = particles[0].z;
  double Vz0f = -1.0*sqrt(pow(Vz0, 2) + 2*g*z0);
  t_base = (Vz0-Vz0f)/g; //Tiempo que tarda una partícula en colisionar con la base

  //Frontera externa:
  
  double t_bext = 0.0, VzN = particles[N-1].Vz, zN = particles[N-1].z;
  
  if (VzN >= 0) //Se calcula el tiempo tbext si la partícula se mueve hacia la frontera exterior 
    {
      double dis = abs(Rf-zN); //Aunque en teoría debería ser siempre positivo, aún falta cuadrar los rangos de las posiciones
      double VzNf = sqrt(pow(VzN, 2)-2*g*dis); //pues dis>0 y es una desaceleración por eso el menos
      t_bext = (VzN-VzNf)/g;
    }
  else{
    t_bext = infty;
  }

  //Asigna los tiempos mínimos de colisión de la primera y ultima partícula a la estructura collisions.
  collisions.t_intb = t_base;
  collisions.t_extb = t_bext;

  //Guarda el tiempo mínimo de colision con las fronteras en el vector times que almacena todos los tiempos de colisiones. Posteriormente en la funcion main se determina el tiempo minimo incluyendo los de las colisiones con las fronteras y las de las partículas vecinas.
  
  times[N-1] = collisions.t_intb;
  times[N] = collisions.t_extb; 
}

//Revisar notación

double Collide_particles_and_evolve_system (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  //Son tres posibles casos: 1) El tiempo mínimo pertenece a una colisión entre partículas vecinas, 2) el tiempo mínimo es una colisión entre la base y la partícula 1 o 3) el tiempo mínimo corresponde a que la última partícula ha sobrepasado la frontera exterior.
  //genera numeros distribuidos uniformemente (Toca editar porque los perfiles de velocidad tienen una distribución maxwelliana)

  vector<double> smallest_time (3, 0.0);
  smallest_time[0] = collisions.t_min; //Tiempo mínimo de las colisiones entre partículas vecinas
  smallest_time[1] = collisions.t_intb; //Tiempo en que la partícula 0 choca con la base
  smallest_time[2] = collisions.t_extb; //Tiempo que la última partícula sale de la frontera exterior
  
  double tt = *min_element(smallest_time.begin(), smallest_time.end());

  movement_equations (particles, tt); // Actualiza el sistema en el tiempo tt mínimo. (Es mejor separar otra función para evolucionar el sistema!!!)
  
  //Caso 1: Revisar muy bien como funciona el cambio de velocidades y phi y theta 
  /*
  if (tt == collisions.t_min)
    {
      std::mt19937 genP(seed);
      std::uniform_real_distribution<double> dis(0, 2*M_PI);
  
      double phi = 2*M_PI*ran.r(), P = dis(genP); //Número aleatorio entre 0 a 2pi
      double theta = acos(sqrt(P));
      
      double Vx1 = ball1.Vx, Vy1 = ball1.Vy, Vz1 = ball1.Vz;
      double Vx2 = ball2.Vx, Vy2 = ball2.Vy, Vz2 = ball2.Vz;
      
      double Ux = Vx2-Vx1, Uy = Vy2-Vy1, Uz = Vz2-Vz1;
      double  U = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
      
      double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);
      
      double Vx1p = (Vx1+Vx2+U*nx)*0.5, Vx2p = (Vx1+Vx2-U*nx)*0.5;
      double Vy1p = (Vy1+Vy2+U*ny)*0.5, Vy2p = (Vy1+Vy2-U*ny)*0.5;
      double Vz1p = (Vz1+Vz2+U*nz)*0.5, Vz2p = (Vz1+Vz2-U*nz)*0.5;
    }
  */
  //Caso 2:

  if (tt == collisions.t_intb)
    {
      particles[0].z = R0; //Es reenviada al sistema desde la base
      particles[0].Vz *= -1.0; //En principio se cambia la velocidad en z, pero la idea es que sea Maxwelliana
    }

  //Caso 3: Se debe escoger si se reenvía desde la base o desde la frontera exterior dependiendo del campo E
  /*
  if(tt == collisions.t_extb)
    {
      
    }
  */

  
  //ball1.Vx=Vx1p; ball1.Vy=Vy1p; ball1.Vz=Vz1p;
  //ball2.Vx=Vx2p; ball2.Vy=Vy2p; ball2.Vz=Vz2p;
  //Integrar...
  //ball1.z+=ball1.Vz*dt; ball2.z+=ball2.Vz*dt;

  return tt;
}

