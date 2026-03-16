#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include "somme.h"

using namespace std;

int main(int argc, char ** argv)
{
  int i,N;
  N = 100000;
  MPI_Init(&argc,&argv);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  double temps_debut, temps_fin;

  int taille;

  cout << "rank = " << rank << "  size = " << size << endl;
  if ( argc > 1 )
    {
      taille = atoi(argv[1]);
      if ( rank == 0 )
	cout << "Taille : " << taille << endl;
    }
  else
    {
      if ( rank == 0 )
	cout << "Il manque la taille en argument" << endl;
      exit(1);
    }

  // Définition de l'objet UneSomme de type Somme par appel au constructeur
  Somme UneSomme(taille); 

  UneSomme.Partitionnement();
  UneSomme.Sommation();

  MPI_Barrier(MPI_COMM_WORLD);
  
  temps_debut = MPI_Wtime();
  for (i=0; i<N; i++)
    {
      UneSomme.Communication_pire();
    }
  temps_fin = MPI_Wtime();
  if ( rank == 0 )
    cout << "Temps de communication pour la plus mauvaise methode : "
	 << temps_fin-temps_debut << " s" << endl;

  UneSomme.Affichage();

  MPI_Barrier(MPI_COMM_WORLD);
  
  
  temps_debut = MPI_Wtime();
  for (i=0; i<N; i++)
    {
      UneSomme.Communication_intermediaire();
    }
  temps_fin = MPI_Wtime();
  if ( rank == 0 )
    cout << "Temps de communication pour la methode intermediaire : "
	 << temps_fin-temps_debut << " s" << endl;
  UneSomme.Affichage();

  MPI_Barrier(MPI_COMM_WORLD);
  
    
  /*
  temps_debut = MPI_Wtime();
  for (i=0; i<N; i++)
    {
      UneSomme.Communication_optimale();
    }
  temps_fin = MPI_Wtime();
  if ( rank == 0 )
    cout << "Temps de communication pour la methode optimale : "
	 << temps_fin-temps_debut << " s" << endl;
  UneSomme.Affichage();

  MPI_Barrier(MPI_COMM_WORLD);
  */
  
  
  temps_debut = MPI_Wtime();
  for (i=0; i<N; i++)
    {
      UneSomme.Communication_Reduce();
    }
  temps_fin = MPI_Wtime();
  if ( rank == 0 )
    cout << "Temps de communication pour la methode MPI_Allreduce : "
	 << temps_fin-temps_debut << " s" << endl;
  UneSomme.Affichage();

  MPI_Barrier(MPI_COMM_WORLD);
 

  MPI_Finalize();
}
