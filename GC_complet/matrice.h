#ifndef MATRICE_H
#define MATRICE_H

/* Fichier contenant les structures de données utilisées pour stocker les matrices :
 * 
 * - matrice_DF : on considère que toutes les lignes de la matrice sont identiques. On ne stocke alors que les n valeurs non nulles de la ligne
 * 
*/

#include <cstdlib>

struct matrice_DF {
    matrice_DF(int lnb_termes_nn, const double* lbuffer, const int* lcolonnes)
    {
        if  ( lnb_termes_nn == 0 )
            exit(EXIT_FAILURE);
        
        nb_termes_nn = lnb_termes_nn; 
        buffer = new double[nb_termes_nn];
        colonnes = new int[nb_termes_nn];
        
        for (int i=0; i<nb_termes_nn; i++)
        {
            buffer[i] = lbuffer[i];
            colonnes[i] = lcolonnes[i];
        }
    }
    
    ~matrice_DF() { delete[] buffer; delete[]  colonnes; } // destructeur
    
    // Opérateurs d'accès : A[i]  renvoie la valeur buffer[i], A(i) renvoie colonnes[i]
    double &operator[](unsigned short int i) const { return buffer[i]; }
    int &operator()(unsigned short int i) const { return colonnes[i]; }
    
    // Les attributs (= les variables de matrice_DF)
    int nb_termes_nn;
    double* buffer;
    int* colonnes;
};

#endif
