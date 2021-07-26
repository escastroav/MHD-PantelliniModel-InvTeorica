#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"
#include "pantellini.h"

TEST_CASE( "Se verifica que las partículas que colisionan primero tengan la misma altura z", "[collide_particles]" ) {

  vector<particle> particles (N);
  collision time_collisions;
  vector<double> times (N+1, 0.0); //En el espacio [N] y [N-1] (último dos espacios) guardará los tiempos de colisión con las fronteras int y ext.
  
  initial_values_particles (particles);

  for (int ii=0; ii<100; ii++) //La primera no se tiene en cuenta 
    {
      sort_particles (particles);
      Get_Collision_Time(particles, times, time_collisions);
      collision_time_boundary (particles, times, time_collisions);
      Collide_particles (particles, times, time_collisions);
      int I = time_collisions.I;
      if ((time_collisions.I != 0) && (time_collisions.I != N+1))
	{
	  REQUIRE(particles[I-1].z == particles[I-2].z);  //Pues I se guar en I-1
	}
      if (time_collisions.I == 0)
	{
	  REQUIRE(particles[0].z == R0);
	}
      
      if (time_collisions.I == N+1)
	{
	  REQUIRE(particles[N].z == Rf);
	}
    }
}

TEST_CASE( "Se verifica que las partículas siempre se encuentren en el dominio", "[evolve_system]" )
{
  vector<particle> particles (N);
  collision time_collisions;
  vector<double> times (N+1, 0.0); //En el espacio [N] y [N-1] (último dos espacios) guardará los tiempos de colisión con las fronteras int y ext.
  
  initial_values_particles (particles);

  for (int ii=0; ii<100; ii++) //La primera no se tiene en cuenta 
    {
      sort_particles (particles);
      Get_Collision_Time(particles, times, time_collisions);
      collision_time_boundary (particles, times, time_collisions);
      Collide_particles (particles, times, time_collisions);
      for (int jj=0; jj<N; jj++)
	{
	  REQUIRE(particles[ii].z >= R0);
	  REQUIRE(particles[ii].z <= Rf);
	}  
    }
}

TEST_CASE( "Se mira que al sumar sobre todo es espacio la densidad de particulas de el número total de partículas N", "[measure_density]" )
{
  vector<particle> particles (N);
  collision time_collisions;
  vector<double> times (N+1, 0.0); //En el espacio [N] y [N-1] (último dos espacios) guardará los tiempos de colisión con las fronteras int y ext.
  
  initial_values_particles (particles);

  for (int ii=0; ii<100; ii++) //La primera no se tiene en cuenta 
    {
      sort_particles (particles);
      Get_Collision_Time(particles, times, time_collisions);
      collision_time_boundary (particles, times, time_collisions);
      Collide_particles (particles, times, time_collisions); 
    }
  int sum = density_of_particles(particles);
  REQUIRE(sum == N);
}
