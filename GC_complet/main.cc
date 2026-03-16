#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "mpi.h"

#include "matrice.h"

using namespace std;

/* *** Résolution d'un système linéaire AX = b en parallèle,
 * - par une méthode directe : décomposition LU de A
 * - par une méthode itérative  : gradient conjugué et gradient conjugué préconditionné
 * 
 * A et b sont issus de la discrétisation par différences finies de l'équation :
 * d^2u/dx^2 + d^2u/dy^2 = -f, soit Laplacien(u) = -f, ou u : [0,1]x[0,1] -> R
 * Le vecteur X contient alors la valeur approchée de u(x,y) en chaque noeud de la grille.
 *
 * En entrée : nombre de points selon x, nombre de points selon y et choix le la méthode (1 GC, 2 GCP)
 * 
*/ 

extern void GC(const matrice_DF& A, const double* b, double* X, double epsilon, int itmax);
extern void GCPC(const matrice_DF& A, const double* b, double* X, double epsilon, int itmax);

int nx=0, ny=0;
double dx=0, dy=0;

int* partitionnement=0; // tableau  d'entiers
int rang=0, nb_proc=0;

void get_coordinates(int indice, double coord[2])
{
     int i_y = indice/nx; // partie entière
    int i_x = indice%nx; // reste de la division : i = i_y*nx + i_x
    
    coord[0] = i_x*dx;
    coord[1] = i_y*dy;
}

void init_scd_membre(int bornes[2], double* b)
{    
    // Ex : b = x*y*(1-x)*(1-y)
    int ind = 0;
    double coord[2];
    for (int i=bornes[0]; i<bornes[1]; i++) // boucle sur i
    {
        get_coordinates(i,coord);
        
        b[ind] = 100*coord[0]*coord[1]*(1-coord[0])*(1-coord[1] );  
        ind++; // Incrémente de 1 <=> ind = ind+1        
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    // Étape 1 : on récupère les données de l'environnement MPI : nombre total de processus et rang du processus
    MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rang);  // nb_proc et rang sont déclarés dans matrice.h
    
    
    // Étape 2 : on récupère les données de l'utilisateur
    if ( ( argc < 4 ) && ( rang == 0) ) // && est l'opérateur logique AND
    {
        printf("Vous devez passer 3 arguments en ligne de commande : \n - Le nombre de points selon x \n - Le nombre de points selon y \n - Le numéro de la méthode de résolution (1 GC, 2 GCP)\n");
        exit(EXIT_FAILURE);
    }
    
    int methode;
    nx = atoi(argv[1]); // nx, ny, dx, dy sont déclarés dans matrice.h
    ny = atoi(argv[2]);
    dx = 1.0/(nx-1); dy = 1.0/(ny-1);
    methode = atoi(argv[3]);
    
    if ( rang == 0 )
        printf("Données lues : nx = %d, ny = %d, méthode = %d\n",nx,ny,methode);
    
    // Étape 3 : partitionnement de la grille de calcul : noeud n0 = (0,0), n1 = (dx,0), ..., n_(nx*ny-1) = (1,1). 
    // Construction du tableau partitionnement de dimension 2*nb_proc. Le noeud i appartient au rang r si partitionnement[2*r] <= i < partitionnement[2*r+1]
    int partitionnement_local[2] = {rang*nx*ny/nb_proc,(rang+1)*nx*ny/nb_proc};
    partitionnement = new int[2*nb_proc]; // Allocation dynamique du tableau partitionnement (2*nb_proc entiers) déclaré dans matrice.h
    MPI_Allgather(partitionnement_local,2,MPI_INT,partitionnement,2,MPI_INT,MPI_COMM_WORLD);
    
    // Étape 4  : déclaration d'une structure matrice adaptée aux différences finies. On suppose que les termes non nuls sont identiques de ligne 
    // à ligne. On doit indique le nombre de termes non nuls par ligne, la valeur de ces termes, et leur numéro de colonne par rapport à la diagonale.
    // Exemple : 0 = terme diagonale, -1 = une colonne avant celle de la diagonale, +1 = une colonne après celle de la diagonale.
    
    double lambda = 1; // conductivité
    double aux_x = lambda/(dx*dx), aux_y = lambda/(dy*dy);
    double lbuffer[] = {2*aux_x+2*aux_y,-aux_x,-aux_x,-aux_y,-aux_y};
    int lcolonnes[] = {0,-1,+1,-nx,+nx}; // Les conditions aux bords seront traitées lors de la résolution (par simplicité)
    matrice_DF A(5,lbuffer,lcolonnes);
    
    // Étape 5 : déclaration du second b membre et du vecteur solution X
    int taille_locale = partitionnement_local[1] - partitionnement_local[0];
    double* b = new double[taille_locale];
    double* X = new double[taille_locale];
    init_scd_membre(partitionnement_local,b);
    for (int i=0; i<taille_locale; i++)
        X[i] = 0;
    
    // Résolution par Gradient Conjugué
    if ( methode == 1 )
        GC(A,b,X,1.0e-4,3000);
    else
        GCPC(A,b,X,1.0e-4,3000);
    
    // Étape finale : écriture d'un fichier au format VTK
    // On rassemble d'abord le champ solution sur le noeud 0
    double* solution  = 0;
    int* deplacements = 0;
    int* nb_data_par_proc = 0;
    
    if ( rang == 0 )
    {
        solution = new double[nx*ny];
        
        deplacements = new int[nb_proc];
        nb_data_par_proc = new int[nb_proc];
        
        for (int i=0; i<nb_proc; i++)
        {
            nb_data_par_proc[i] = partitionnement[2*i+1] - partitionnement[2*i];
        }
        
        deplacements[0] = 0;
        for (int i=1; i<nb_proc; i++)
        {
            deplacements[i] = deplacements[i-1] + nb_data_par_proc[i-1];
        }
    }
    
    int nb_data = partitionnement_local[1] - partitionnement_local[0];
    MPI_Gatherv(X,nb_data,MPI_DOUBLE,solution,nb_data_par_proc,deplacements,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    delete[] nb_data_par_proc;
    delete[] deplacements;
    
    if ( rang == 0 )
	{    
        printf("Écriture du fichier vtk");
        stringstream nom;

        nom << "TP_CHP";
        nom << ".vtk" << '\0';
        ofstream fic(nom.str().c_str());

        fic << "# vtk DataFile Version 2.0" << endl;
        fic << "Résolution d'un système linéaire" << endl;
        fic << "ASCII" << endl;
        fic << "DATASET STRUCTURED_POINTS" << endl;
        fic << "DIMENSIONS " << nx << "  " << ny << "  1 " << endl;
        fic << "ORIGIN 0 0 0" << endl;
        fic << "SPACING " << dx << "  " << dy << "  1" << endl;
        fic << "POINT_DATA " << nx*ny << endl;
        fic << "SCALARS solution float" << endl;
        fic << "LOOKUP_TABLE default" << endl;
    
        for (int i=0; i<nx*ny; i++)
        {
            fic << solution[i] << endl;
        }
        fic.close();
    }
    
    delete[] solution;
      
    delete [] X;
    delete [] b;
    delete [] partitionnement;
    MPI_Finalize();

   return EXIT_SUCCESS;
}
