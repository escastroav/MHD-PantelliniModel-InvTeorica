#include <iostream>
#include <cmath>
//#include "Random64.h"
#include "pantellini.h"
#include <iomanip>

using namespace std;

int main()
{
  cout << fixed;cout.precision(20);
  Crandom rand64(1);
  Ball balls[N];
  Collider colls;

  double tColl = 0, time = 0; 
  int indexColl = 0;
  int collisions = 0;

  bool areBallsSorted = true;

  //1. Inicializar las posiciones de las particulas.
  colls.InitPositions(balls, rand64);
  //InicieAnimacion();
  //cout << "time" << "\t" << "I" << "\t" << "tmin" << "\t" << "z_I" << "\t" << "Vz_I" << "\n" ;
  while(time <= 1000*tf)
    {    
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
	colls.CollideBalls(balls[indexColl], balls[indexColl-1], rand64, tColl);
	//balls[indexColl-1].Integrate(tColl);
	//balls[indexColl].Integrate(tColl);
      }
      //6. Imprimir y actualizar el siguiente paso.
      /*InicieCuadro();
      for(int i=0;i<N;i++)
	  balls[i].PrintBall();
	  TermineCuadro();*/
      
	  time += tColl;
	  collisions++;
	  //cout << time << "\t" << indexColl << "\t" << tColl << "\t" << balls[indexColl].GetZ() << "\t" << balls[indexColl].GetVz() << "\n" ;
    }
  colls.MeasureDensity(balls);
  //for(int i=0;i<N;i++)cout << balls[i].GetI() << "\t" << balls[i].GetZ() << "\n";
  //cout << collisions << endl;
  return 0;
}
