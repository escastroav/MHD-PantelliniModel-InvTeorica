#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <chrono>
#define _USE_MATH_DEFINES

using namespace std;
//Definición de constantes

const double size = 1;
const int N = 50;//*size;
const double R0 = 5.0;
const double Rf = 50*R0;
const int seed = 0;
const double h = N*R0;
const double g = 3*R0;
const double H = 0.2; // kT/mg=rmax/5
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
double Collide_particles (vector<particle> & particles, vector<double> & times, collision & collisions);
void Maxwellian_distribution (vector<particle> & particles);
void evolve_system (vector<particle> & particles, vector<double> & times, collision & collisions);
void print_some_stuff (vector<particle> & particles, vector<double> & times, collision & collisions);
void density_of_particles (vector<particle> particles);

int main ()
{
  std::cout.precision(6);
  std::cout.setf(std::ios::scientific);

  //Inicialización de todos los objetos necesarios
  vector<particle> particles (N);
  collision time_collisions;
  vector<double> times (N+1, 0.0); //En el espacio [N] y [N-1] (último dos espacios) guardará los tiempos de colisión con las fronteras int y ext.

  print_some_stuff (particles, times, time_collisions);

  return 0;
}

void initial_values_particles (vector<particle> & particles) //Tiene como objetivo dar valores iniciales aleatorios a las partículas
{
  //Enumera las partículas
  for (int ii=0; ii<N; ii++)
    {
      particles[ii].id = ii+1;
      //genera numeros distribuidos uniformemente (Toca editar porque los perfiles de velocidad tienen una distribución maxwelliana)
      std::mt19937 gen(ii);
      std::uniform_real_distribution<double> dis(0, 25.0); //Rango de velocidades arbitrario
      std::uniform_real_distribution<double> disp(R0, Rf); //Rangos del dominio de la simulación
      particles[ii].z = disp(gen);
      particles[ii].Vx = dis(gen);
      particles[ii].Vy = dis(gen);
      particles[ii].Vz = 0; //dis(gen);
    }
  sort_particles (particles);
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
      particles[ii].id = ii+1;
    }
}

void Get_Collision_Time(vector<particle> & particles, vector<double> & times, collision & collisions)
{
  
  for (int ii=0; ii<(N-1); ii++) //Pues el N-1 asegura que no se guarde ningun tiempo en el último espacio (designado para tiempos con las fronteras)
    {
      double vrel = (particles[ii+1].Vz-particles[ii].Vz);
      if((abs(vrel) > eps) && (particles[ii+1].z != particles[ii].z)) //Las particulas no se encuentran en la misma posición
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
	  collisions.I = ii+2;  //Pues el elemento times[ii] te da tiempo entre ii y ii+1, significa que I=ii+2 y I-1=ii+1 pues ii empieza en cero y no en 1. (Suponga que el tmin es entre la particula 1 y 2, I deberia ser 2 y no 1 pues en nuestra notacion no existe la particula cero, pero en los vectores el cero le corresponde a la particula 1)
	  break;
	}
    }
}

void collision_time_boundary (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  //Frontera interna
  
  double t_base = infty, Vz0 = particles[0].Vz, z0 = particles[0].z;
  
  if (z0 != R0){  //La partícula no se encuentre en la base pues daría un tiempo igual a cero y eso significa que ya colisionó antes
    double Vz0f = -1.0*sqrt(pow(Vz0, 2) + 2*g*abs(z0-R0)); //En teoría z0 siempre debería ser mayor a R0.
    t_base = (Vz0-Vz0f)/g; //Tiempo que tarda una partícula en colisionar con la base
  }
  
  //Frontera externa:
  
  double t_bext = infty, VzN = particles[N-1].Vz, zN = particles[N-1].z;
  
  if (zN != Rf){ //Verifica que la partícula no se encuentre en la frontera exterior
    double dis = abs(Rf-zN); //Aunque en teoría debería ser siempre positivo, aún falta cuadrar los rangos de las posiciones
    double deter = pow(VzN, 2)-2*g*dis; //Término dentro de la raíz el cual debe ser positivo si la partícula puede alcanzar la fron. ext
    
    if ((VzN >= 0) && (deter >= 0)) //Se calcula el tiempo tbext si la partícula se mueve hacia la frontera exterior 
      {
	double VzNf = sqrt(deter); //pues dis>0 y es una desaceleración por eso el menos
	t_bext = (VzN-VzNf)/g;
      }
    else{
      t_bext = infty;
    }
  }
  
  //Asigna los tiempos mínimos de colisión de la primera y ultima partícula a la estructura collisions.
  collisions.t_intb = t_base;
  collisions.t_extb = t_bext;

  //Guarda el tiempo mínimo de colision con las fronteras en el vector times que almacena todos los tiempos de colisiones. Posteriormente en la funcion main se determina el tiempo minimo incluyendo los de las colisiones con las fronteras y las de las partículas vecinas.

  times[N-1] = collisions.t_intb;
  times[N] = collisions.t_extb; 
}

double Collide_particles (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  //Son tres posibles casos: 1) El tiempo mínimo pertenece a una colisión entre partículas vecinas, 2) el tiempo mínimo es una colisión entre la base y la partícula 1 o 3) el tiempo mínimo corresponde a que la última partícula ha sobrepasado la frontera exterior.
  //genera numeros distribuidos uniformemente (Toca editar porque los perfiles de velocidad tienen una distribución maxwelliana)

  vector<double> smallest_time (3, 0.0);
  
  smallest_time[0] = collisions.tmin; //Tiempo mínimo de las colisiones entre partículas vecinas
  smallest_time[1] = collisions.t_intb; //Tiempo en que la partícula 0 choca con la base
  smallest_time[2] = collisions.t_extb; //Tiempo que la última partícula sale de la frontera exterior

  double tt = *min_element(smallest_time.begin(), smallest_time.end());

  movement_equations (particles, tt); // Actualiza el sistema en el tiempo tt mínimo.
  
  //Caso 1: Revisar muy bien como funciona el cambio de velocidades y phi y theta 

  if (tt == collisions.tmin)
    {
      int I = collisions.I;
      
      std::mt19937 genP(seed); //Recomendable que la semilla no sea la misma para cuando es sistema tenga que evolucionar muchas veces.
      std::uniform_real_distribution<double> disP(0, 1.0);
      std::mt19937 genphi(seed);
      std::uniform_real_distribution<double> disphi(0, 2*M_PI);
      
      double phi = disphi(genphi), P = disP(genP); //Número aleatorio entre 0 a 2pi para phi y 0 a 1 para P
      double theta = acos(sqrt(P));

      double Vx1 = particles[I-1].Vx, Vy1 = particles[I-1].Vy, Vz1 = particles[I-1].Vz; 
      double Vx2 = particles[I-2].Vx, Vy2 = particles[I-2].Vy, Vz2 = particles[I-2].Vz; 
      
      double Ux = Vx2-Vx1, Uy = Vy2-Vy1, Uz = Vz2-Vz1;
      double U = sqrt(Ux*Ux+Uy*Uy+Uz*Uz);
      
      double nx = cos(phi)*sin(theta), ny = sin(phi)*sin(theta), nz = cos(theta);
      //¿Cómo se calculan estas velocidades? Revisar!!!!!!!!!!!
      double Vx1p = (Vx1+Vx2+U*nx)*0.5, Vx2p = (Vx1+Vx2-U*nx)*0.5;
      double Vy1p = (Vy1+Vy2+U*ny)*0.5, Vy2p = (Vy1+Vy2-U*ny)*0.5;
      double Vz1p = (Vz1+Vz2+U*nz)*0.5, Vz2p = (Vz1+Vz2-U*nz)*0.5;
      //Se actualizan las velocidades de la Partícula I después de la colisión
      particles[I-1].Vx = Vx1p;
      particles[I-1].Vy = Vy1p;
      particles[I-1].Vz = Vz1p;
      //Se actualizan las velocidades de la Partícula I-1 después de la colisión
      particles[I-2].Vx = Vx2p;
      particles[I-2].Vy = Vy2p;
      particles[I-2].Vz = Vz2p;
      //Se actualizan las posiciones de las partículas I e I-1. es importante definirlos así
      particles[I-1].z = particles[I-2].z;
    }
  
  //Caso 2:

  if (tt == collisions.t_intb)
    {
      particles[0].z = R0; //Es reenviada al sistema desde la base
      particles[0].Vz *= -1.0; //En principio se cambia la velocidad en z, pero la idea es que sea Maxwelliana
    }

  //Caso 3: Se debe escoger si se reenvía desde la base o desde la frontera exterior dependiendo del campo E
  
  if(tt == collisions.t_extb)
    {
      particles[N-1].z = Rf; //Es reenviada al sistema desde la frontera exterior.
      particles[N-1].Vz *= -1.0; //En principio se cambia la velocidad en z, pero la idea es que sea Maxwelliana
    }
  
  return tt;
}

void evolve_system (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  initial_values_particles (particles);
  /*
  std::cout<<"Se imprimen para las 5 partículas sus valores para id, z, Vx, Vy, Vz después de darles unos valores iniciales:"<<std::endl;
  
  for (int jj=0; jj<N; jj++)
    {
      std::cout<<particles[jj].id<<"\t"<<particles[jj].z<<"\t"<<particles[jj].Vx<<"\t"<<particles[jj].Vy<<"\t"<<particles[jj].Vz<<std::endl;
    }
  */
  for (int ii=0; ii<20; ii++) //Hacemos un paso de evolución más
    {
      sort_particles (particles);
      Get_Collision_Time(particles, times, collisions);
      collision_time_boundary (particles, times, collisions);
      /*
      std::cout<<"Se imprimen los tiempos mínimos de la colisión "<<ii+1<<" (entre partículas y fronteras)"<<std::endl;
      std::cout<<collisions.I<<"\t"<<collisions.tmin<<"\t"<<collisions.t_intb<<"\t"<<collisions.t_extb<<"\t"<<*min_element(times.begin(), times.end())<<std::endl;
      */
      Collide_particles (particles, times, collisions);
      /*
      std::cout<<"Se imprimen para las 5 partículas sus valores para id, z, Vx, Vy, Vz después de evolucionar por "<<ii+1<<" vez el sistema"<<std::endl;

      for (int jj=0; jj<N; jj++)
      {
      std::cout<<particles[jj].id<<"\t"<<particles[jj].z<<"\t"<<particles[jj].Vx<<"\t"<<particles[jj].Vy<<"\t"<<particles[jj].Vz<<std::endl;
      }
      */
    }
  
  std::cout<<"Se imprimen para las 5 partículas sus valores para id, z, Vx, Vy, Vz después de darles unos valores iniciales:"<<std::endl;
  for (int jj=0; jj<N; jj++)
    {
      std::cout<<particles[jj].id<<"\t"<<particles[jj].z<<"\t"<<particles[jj].Vx<<"\t"<<particles[jj].Vy<<"\t"<<particles[jj].Vz<<std::endl;
    }
}


void print_some_stuff (vector<particle> & particles, vector<double> & times, collision & collisions)
{
  evolve_system (particles, times, collisions);
  density_of_particles (particles);
}

void density_of_particles (vector<particle> particles)
{
  double dz = 5.0 , deltaZ = Rf-R0; //deltaZ = 245
  int M = deltaZ/dz; //M=49  
  std::vector<int> density (M, 0.0);
  int sum_n = 0;

  cout << "Altura z promedio"<< "\t" << "Densidad de partículas" << "\n";
  
  for(int ll=0; ll<M; ll++)
    {
      double zi = R0;
      
      zi += ll*dz;
      
      for(int ii = 0; ii<N; ii++)
	{	  
	  if((particles[ii].z > zi) && (particles[ii].z < zi+dz))
	    density[ll]++;
	}
      
      if(density[ll] != 0)
	{
	  cout << zi+dz/2 << "\t" << density[ll] << "\n"; //Escoge a zi el valor medio de altura del intervalo escogido
	  sum_n += density[ll];
	}
    }
  
  cout << "Número de partículas N resultado de sumar por cada intervalo de z las densidades: "<<sum_n<<"\n";
}

/*
Posibles problemas a tener en cuenta:
- Si una partícula está en la base o en la frontera exterior, no la tenemos en cuenta para el paso de evolución, pero y sus partículas vecinas? puede darse el caso en que la segunda esté mas cerca de la base que cualesquiera dos partículas colisionen y lo mismo para la otra frontera. Si solo tenemos la coordenada z no debería haber inconveniente, pues primero chocan las primeras dos antes de que la segunda toca la base.

-Falta Hacer la función de densidad de partículas.

*/

/*
void Maxwellian_distribution (vector<particle> & particles)
{
  unsigned seed2 = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);
  
  // Boltzmann factor times temperature
  const double k_T = 0.1;
  
  // setup the Maxwell distribution, i.e. gamma distribution with alpha = 3/2
  gamma_distribution<double> maxwell(3./2., k_T);
  
  // generate Maxwell-distributed values
  for (int ii = 0; ii < 10000; ii++) {
    particles[ii].Vz = maxwell(generator);
    cout <<maxwell(generator) << endl;
  }
}
*/
