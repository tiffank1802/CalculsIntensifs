/* Méthode du gradient conjugué.

*/

#include <cstdio>
#include <cmath>

#include "matrice.h"
#include "mpi.h"

extern int nx, ny, dx, dy;
extern int*  partitionnement; // tableau d'entiers
extern int rang, nb_proc;

extern void get_coordinates(int indice, double coord[2]); // définie dans main.cc

void ProdMatVecSeq(const matrice_DF& A, const double* X, double* produit, int taille)
{
    double coord[2];
    double precision = 1.0e-8;
    for (int i=0; i<taille; i++)  // boucle sur les lignes de A et produit
    {
        get_coordinates(i,coord);
        produit[i] = 0;
        for (int j=0; j<A.nb_termes_nn; j++)
        {
            if ( ( coord[0] > precision ) && ( coord[0] < 1-precision ) && ( coord[1] > precision ) && ( coord[1] < 1-precision ) ) // condition au bord homogène
            {
                int col = i + A(j);
                if ( ( col >= 0 ) && ( col < taille ) )  // condition satisfaite puisque l'on n'est pas en y=0 ou y=1
                    produit[i] = produit[i] + A[j]*X[col];
            }            
        }      
    }    
}

void ProdMatVecPara(const matrice_DF& A, const double* X, double* produit, int taille)
{
    // Étape 1 : le processus Pn envoie au processus Pn+1 ses nx dernières valeurs de X. Pn+1 les reçoit.
    // Étape 2 : Pn+1 envoie à Pn ses nx premières valeurs de X. Pn les reçoit.
    
    double* X_de_Pn = new double[nx];
    double* X_de_Pnp1 = new double[nx];
    
    // 1.
    
    if ( taille - nx < 0 )
    {
        printf("Processus %d : pas assez de données sur ce processus\n",rang);
        exit(EXIT_FAILURE);
    }
    
    int tag1 = 0;
    if ( rang < nb_proc-1 )
    {
        const double* buffer = X + taille - nx;
        MPI_Send(buffer,nx,MPI_DOUBLE,rang+1,tag1,MPI_COMM_WORLD);        
    }
    
    if ( rang > 0 )
    {
        MPI_Recv(X_de_Pn,nx,MPI_DOUBLE,rang-1,tag1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);        
    }
    
    // 2.
    
    int tag2 = 1;
    if ( rang > 0 )
    {
        MPI_Send(X,nx,MPI_DOUBLE,rang-1,tag2,MPI_COMM_WORLD);
    }
    
    if ( rang < nb_proc-1 )
    {
        MPI_Recv(X_de_Pnp1,nx,MPI_DOUBLE,rang+1,tag2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);        
    }
    
    int bornes[2] = {partitionnement[2*rang],partitionnement[2*rang+1]};        
    double coord[2];
    double precision = 1.0e-8;
    for (int i=0; i<taille; i++)  // boucle sur les lignes de A et produit
    {
        get_coordinates(i+bornes[0],coord);
        produit[i] = 0;
        for (int j=0; j<A.nb_termes_nn; j++)
        {
            if ( ( coord[0] > precision ) && ( coord[0] < 1-precision ) && ( coord[1] > precision ) && ( coord[1] < 1-precision ) ) // condition au bord homogène
            {
                int col = i + A(j);
                double val_X;
                if ( col < 0 )
                {
                    val_X = X_de_Pn[col+nx];
                }
                else
                {
                    if ( col >= taille )
                    {
                         val_X = X_de_Pnp1[col-taille];
                    }
                    else
                        val_X = X[col];
                }
                
                produit[i] = produit[i] + A[j]*val_X;
            }            
        }      
    } 
      
    delete [] X_de_Pn;
    delete [] X_de_Pnp1;
}

double ProdScalaire(const double* vec1, const double* vec2, int taille)
{
    double produit_local = 0;
    for (int i=0; i<taille; i++)
        produit_local = produit_local + vec1[i]*vec2[i];
    
    double produit = 0;
    MPI_Allreduce(&produit_local,&produit,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    return produit;
}

void GC(const matrice_DF& A, const double* b, double* X, double epsilon, int itmax)
{
    double t0 = MPI_Wtime();
    int taille_locale = partitionnement[2*rang+1] - partitionnement[2*rang];
       
    // Algorithme du gradient conjugué    
    double* ri = new double[taille_locale];
    double* di = new double[taille_locale];
    double* produit = new double[taille_locale];
    
    ProdMatVecPara(A,X,produit,taille_locale);
    for (int k=0; k<taille_locale; k++)
    {
        double val = b[k] - produit[k];        
        ri[k] = val;
        di[k] = val;
    }
    
    int it=0;
    
    double ri_ri = ProdScalaire(ri,ri,taille_locale);
    
    if ( rang == 0 )
        printf("Résidu initial : %.3e\n",sqrt(ri_ri));
    
    while ( ( sqrt(ri_ri) > epsilon ) && ( it < itmax ) )  // convergence sur le résidu relatif
    {
        ProdMatVecPara(A,di,produit,taille_locale);
        double alphai = ProdScalaire(ri,di,taille_locale)/ProdScalaire(di,produit,taille_locale);
        for (int k=0; k<taille_locale; k++)
        {
            X[k] = X[k] + alphai*di[k];
            ri[k] = ri[k] - alphai*produit[k];            
        }
        
        double rip1_rip1 = ProdScalaire(ri,ri,taille_locale);
        if ( sqrt(rip1_rip1) > epsilon )
        {
            double betai = rip1_rip1/ri_ri;
            for (int k=0; k<taille_locale; k++)
            {
                di[k] = ri[k] + betai*di[k];
            }            
        }
        
        ri_ri = rip1_rip1;
        
        if ( (it % 10 == 0 ) && (rang == 0) )
            printf("Fin de l'itération %d, résidu %.3e\n",it,sqrt(ri_ri));
        
        it++;
    }
    
    if ( rang == 0 )
        printf("Résidu final %.3e, atteint en %d itérations et %f s.\n",sqrt(ri_ri),it,MPI_Wtime()-t0);
    
    
    delete [] produit;
    delete [] di;
    delete [] ri;
}

void GCPC(const matrice_DF& A, const double* b, double* X, double epsilon, int itmax)
{
    double t0 = MPI_Wtime();
    int taille_locale = partitionnement[2*rang+1] - partitionnement[2*rang];
    
    // Application du préconditionnement
    double diag = A[0]; 
    
    // Algorithme du gradient conjugué    
    double* ri = new double[taille_locale];
    double* zi = new double[taille_locale];
    double* di = new double[taille_locale];
    double* produit = new double[taille_locale];
    
    ProdMatVecPara(A,X,produit,taille_locale);
    for (int k=0; k<taille_locale; k++)
    {
        double val = b[k] - produit[k];        
        ri[k] = val;
        zi[k] = val/diag;
        di[k] = zi[k];
    }
    
    int it=0;
    
    double ri_ri = ProdScalaire(ri,ri,taille_locale);
    double zi_ri = ProdScalaire(zi,ri,taille_locale);
    
    if ( rang == 0 )
        printf("Résidu initial : %.3e\n",sqrt(ri_ri));
    
    while ( ( sqrt(ri_ri) > epsilon ) && ( it < itmax ) ) // convergence sur le résidu relatif
    {
        ProdMatVecPara(A,di,produit,taille_locale);
        double alphai = ProdScalaire(ri,zi,taille_locale)/ProdScalaire(di,produit,taille_locale);
        for (int k=0; k<taille_locale; k++)
        {
            X[k] = X[k] + alphai*di[k];
            ri[k] = ri[k] - alphai*produit[k];            
        }
        
        double rip1_rip1 = ProdScalaire(ri,ri,taille_locale);
        if ( sqrt(rip1_rip1) > epsilon )
        {
            for (int k=0; k<taille_locale; k++)
            {
                zi[k] = ri[k]/diag;
            }
            double aux = zi_ri;
            zi_ri = ProdScalaire(zi,ri,taille_locale);            
            double betai = zi_ri/aux;
            for (int k=0; k<taille_locale; k++)
            {
                di[k] = zi[k] + betai*di[k];
            }            
        }
        
        ri_ri = rip1_rip1;
        
        if ( (it % 10 == 0 ) && (rang == 0) )
        {
            printf("Fin de l'itération %d, résidu : %.3e\n",it,sqrt(ri_ri));
        } 
        
        it++;
    }
    
    if ( rang == 0 )
    {
         printf("Résidu final : %.3e, atteint en %d itérations et %f s.\n",sqrt(ri_ri),it,MPI_Wtime()-t0);
    }    
    
    delete [] produit;
    delete [] di;
    delete [] zi;
    delete [] ri;
}

