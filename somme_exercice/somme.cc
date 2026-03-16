#include <iostream>
#include <iomanip>
#include "somme.h"
#include "mpi.h"

using namespace std;

Somme::Somme(int taille)
{
	taille_totale = taille;
	debut_local = 0;
	fin_local = 0;
	sum = 0;
	sum_p = 0;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
}

void Somme::Partitionnement()
{
	debut_local = taille_totale*rank/nb_proc + 1;
	fin_local = taille_totale*(rank+1)/nb_proc;

	cout << "Proc : " << rank << ", debut = " << debut_local << ", fin = " << fin_local << endl;
}

void Somme::Affichage()
{
	cout << "rang " << rank << ", la somme des entiers de 1 Ã  " << taille_totale << " vaut : " << setprecision(30) << sum << endl;
}


void Somme::Sommation()
{
	// Somme partielle locale au processeur

	int i;
	sum_p = 0;
	for (i=debut_local; i<=fin_local; i++)
	{
		sum_p += i;
	}
}

void Somme::Communication_pire()
{
	// Methode 1 MPI pour faire la somme sur tous les processeurs : pire

	int i;
	for (i=0; i<nb_proc; i++)// pour tous les processus
	{
		if ( i != rank )// pour les processus autres que celui de l'Ãmetteur, il envoie sa somme  partielle aux autres processus
			MPI_Send(&sum_p,1,MPI_INT,i,1,MPI_COMM_WORLD);
	}

	sum = sum_p;
	for (i=0; i<nb_proc; i++)// parcours tous les processsus
	{
		if ( i!= rank ) // pour les processus autres que celui de l'emetteur,il recoit les somme patielles
		{
			int sum_aux;
			MPI_Status status;
			MPI_Recv(&sum_aux,1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
			sum += sum_aux;
		}
	}

}

void Somme::Communication_intermediaire()
{
	// Methode 1 MPI pour faire la somme sur tous les processeurs : l'intermediaire

	if (rank != 0)// si le processus present n'est pas de rang 0, il envoie sa somme partielle au processus de rang 0
	{
		MPI_Send(&sum_p,1,MPI_INT,0,1,MPI_COMM_WORLD); // envoie la somme partielle au processus de rang 0
	}
	else // sinon(s'il est de rang 0, il  reÃoit les sommes partielles des autres processus de rang different de 0 et effectue la somme totale
	{
		sum = sum_p;
		int j;
		for (j=1; j<nb_proc; j++)
		{
			int sum_aux;
			MPI_Status status;
			MPI_Recv(&sum_aux,1,MPI_INT,j,1,MPI_COMM_WORLD,&status);// reÃoit les sommes partielles des autres processus
			sum += sum_aux; // effectue la somme totale
		}

		for (j=1; j<nb_proc; j++)
		{
			MPI_Send(&sum,1,MPI_INT,j,1,MPI_COMM_WORLD); // envoie la somme totale aux autres processus autres que lui-mÃme puisqu'il est de rang 0 (car a effectuÃ© la somme des sommes partielles)
		}
	}

	if (rank != 0) // si le processus actuel n'est pas de rang 0, il reÃoit la somme totale du processus de rang 0
	{
		MPI_Status status;
		MPI_Recv(&sum,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
	}

}

void Somme::Communication_optimale()
{
	// Methode 1 MPI pour faire la somme sur tous les processeurs : l'optimale

	MPI_Status status;

	int sum_aux;

	int inc = 1;
	int dernier = nb_proc - 1;
	while ( inc <= nb_proc/2 )
	{
		if ( ( rank % (2*inc) == inc ) && ( rank <= dernier ) )
		{
			MPI_Send(&sum_p,1,MPI_INT,rank-inc,1,MPI_COMM_WORLD);
		}
		else if ( ( rank + inc <= dernier ) && ( rank % (2*inc) == 0 ) )
		{
			MPI_Recv(&sum_aux,1,MPI_INT,rank+inc,1,MPI_COMM_WORLD,&status);
			sum_p += sum_aux ;
		}

		// Cas du dernier processus, si il est pair

		if ( dernier % (2*inc) == inc )
		{
			dernier -= inc;
		}
		else
		{
			if ( rank == dernier )
			{
				MPI_Send(&sum_p,1,MPI_INT,0,1,MPI_COMM_WORLD);
			}

			if ( rank == 0 )
			{
				MPI_Recv(&sum_aux,1,MPI_INT,dernier,1,MPI_COMM_WORLD,&status);
				sum_p += sum_aux;
			}

			dernier -= 2*inc;
		}

		inc *= 2;
	}

	if ( rank == 0 )
	{
		sum = sum_p;
	}

	// Distribution dans l'autre sens
/*
	inc /= 2;
	dernier = 0;
	while ( inc >= 1 )
	{
		if ( ( rank % (2*inc) == 0 ) && ( rank <= dernier ) && ( rank + inc < nb_proc ))
		{
			...
		}
		else if ( ( rank - inc <= dernier ) && ( rank % (2*inc) == inc ) )
		{
			...
			sum = sum_aux ;
		}

		dernier += inc;
		if ( dernier + inc < nb_proc )
		{
			dernier += inc;

			if ( rank == 0 )
			{
				...
			}

			if ( rank == dernier )
			{
				...
				sum = sum_aux;
			}
		}

		inc /= 2;
	}
*/
}


void Somme::Communication_Reduce()
{
	// Methode "integree" MPI pour faire la somme sur tous les processeurs

	int sum_aux = 0.0;

	// MPI_Reduce(&sum_p,&sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Allreduce(&sum_p,&sum_aux,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	sum = sum_aux;
}
