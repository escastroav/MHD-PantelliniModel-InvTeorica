#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"
#include "Random64.h"

TEST_CASE( "Se revisa que las probabilidades no exceda", "[factorial]" ) {

  
  /*
  unsigned int N=64;
  
  float prob=0.0;
  double * larg_cl= new double[sample_size] {0.0}; //array para los clusters percolantes mas grandes de cada matriz (normalizado)
  double * percol_a= new double[sample_size] {0.0};//array para la percolacion de matrices
  std::vector<statistics> stat_vect;

    while (prob <= 1.0)
    {
      if (prob < 0.55){
	  prob += 0.055;
      }
      if ((prob>=0.55) && (prob<=0.65)){
	prob += 0.01;
      }
      if (prob>0.65){
	prob += 0.035;
	  }
	    if (prob > 1.0){
	break;
      }
	
      for (int seed=0; seed<sample_size; seed++)
	{
	  std::vector<cluster_attributes> cl_att_vect; //vector de atributos del cluster
	  cluster_series_generate(seed, N, prob, cl_att_vect, percol_a);
	  larg_cl[seed]= (largest_perc_cluster(cl_att_vect)/(N*N));
	  std::vector<cluster_attributes>().swap(cl_att_vect);
      
	}
      gen_stat(larg_cl, percol_a, sample_size, prob, stat_vect);

      REQUIRE(gen_stat (double * larg_cl, double * percol_a, int sample_size, 1.0 , std::vector<statistics> & stat_vect) == 1.0);
      REQUIRE(gen_stat (double * larg_cl, double * percol_a, int sample_size, 0.0 , std::vector<statistics> & stat_vect) == 0.0);
    }

  std::vector<statistics>().swap(stat_vect);
  delete [] larg_cl;
  delete [] percol_a;
  return 0;
  */
  
}
