#  Calculs Intensifs — Programmation Parallèle avec MPI en C++

> **Projet académique — École Centrale Lyon ENISE**
> Auteurs : Kevin TONGUE · Gaston T. KAMDEM
> Cluster utilisé : **Centaure** (École des Mines de Saint-Étienne) — 656 cœurs, architecture MIMD

Travaux pratiques de calcul haute performance sur cluster de calcul, utilisant la librairie **MPI** (Message Passing Interface) en **C++**. Le projet couvre la prise en main de MPI, le calcul parallèle sur plusieurs cœurs, et la résolution d'un système linéaire éléments finis par la méthode du **gradient conjugué** parallélisée.

---

## Table des matières

- [Contexte](#contexte)
- [Structure du projet](#structure-du-projet)
- [Prérequis](#prérequis)
- [Compilation & exécution](#compilation--exécution)
- [Description des exercices](#description-des-exercices)
- [Performances mesurées](#performances-mesurées)

---

## Contexte

Les méthodes numériques modernes (éléments finis, dynamique des fluides, etc.) nécessitent des ressources de calcul considérables. Ces TP ont été réalisés sur le cluster **Centaure** de l'École des Mines de Saint-Étienne, qui dispose de **656 cœurs** répartis sur plusieurs nœuds physiques interconnectés par un réseau haute vitesse.

L'architecture utilisée est de type **MIMD** (Multiple Instructions Multiple Data) : chaque processus exécute son propre programme sur ses propres données, communiquant avec les autres via **MPI**.

---

## Structure du projet

```
CalculsIntensifs/
│
├── hello_exercice/          # Prise en main de MPI — Hello World parallèle
│   ├── principal.cc
│   └── Makefile
│
├── somme_exercice/          # Somme des n premiers entiers sur plusieurs cœurs
│   ├── principal.cc         # Programme principal
│   ├── somme.cc             # Implémentation de la classe Somme
│   ├── somme.h              # Déclaration de la classe
│   └── Makefile
│
└── GC_complet/              # Gradient Conjugué parallèle — résolution système linéaire EF
    ├── principal.cc
    ├── *.cc / *.h
    └── Makefile
```

---

## Prérequis

| Outil | Version |
|---|---|
| Compilateur C++ | `g++` ≥ 9 ou `mpicxx` |
| MPI | OpenMPI ≥ 4.x ou MPICH |
| Make | GNU Make |

Installation MPI (Linux) :
```bash
sudo apt install mpich        # ou
sudo apt install libopenmpi-dev openmpi-bin
```

---

## Compilation & exécution

```bash
# Compiler un exercice
cd somme_exercice
make

# Lancer sur 4 cœurs
mpirun -np 4 ./somme

# Lancer sur le cluster (via scheduler)
mpirun -np 16 ./GC_complet/gradient_conjugue
```

---

## Description des exercices

### 1. `hello_exercice` — Fondamentaux MPI

Introduction aux primitives de base de la librairie MPI :

```cpp
#include <mpi.h>

MPI_Init(&argc, &argv);              // Initialisation
MPI_Comm_size(MPI_COMM_WORLD, &p);  // Nombre total de processus
MPI_Comm_rank(MPI_COMM_WORLD, &id); // Rang (identifiant) du processus courant
MPI_Finalize();                      // Fermeture
```

Chaque processus s'identifie et affiche son rang au sein du communicateur `MPI_COMM_WORLD`. Introduit les notions de **communication point-à-point** et **collective**.

---

### 2. `somme_exercice` — Somme parallèle des n premiers entiers

Calcul distribué de `Σ(1..n)` sur `p` cœurs avec **partitionnement des données** : chaque processus calcule une somme partielle sur son intervalle `[debut_local, fin_local]`, puis les résultats sont agrégés.

Trois stratégies de communication comparées :

| Méthode | Principe | Complexité |
|---|---|---|
| **Naïve** | Chaque processus envoie à tous les autres (`MPI_Send` / `MPI_Recv`) | O(p²) |
| **Intermédiaire** | Tous envoient au processus 0 qui agrège | O(p) |
| **Optimale** | Réduction par arbre binaire (pair/impair) | O(log p) |

Fonctions MPI de réduction natives :
```cpp
MPI_Reduce(...)    // Réduit vers le processus 0 uniquement
MPI_Allreduce(...) // Réduit et redistribue à tous les processus
MPI_Barrier(MPI_COMM_WORLD); // Synchronisation entre processus
```

Mesure du temps avec `MPI_Wtime()`.

---

### 3. `GC_complet` — Gradient Conjugué parallèle (méthode de Krylov)

Résolution d'un système linéaire `Ax = b` typique des problèmes **éléments finis**, avec `A` symétrique définie positive. La méthode du gradient conjugué minimise la fonctionnelle quadratique :

```
J(x) = ½(Ax, x) − (b, x)
```

**Algorithme itératif :**

```
Initialisation : x⁽⁰⁾, r⁽⁰⁾ = b − Ax⁽⁰⁾, d⁽⁰⁾ = r⁽⁰⁾
Pour i = 0, 1, 2, ...
    α⁽ⁱ⁾ = (r⁽ⁱ⁾, d⁽ⁱ⁾) / (d⁽ⁱ⁾, Ad⁽ⁱ⁾)   ← pas de descente
    x⁽ⁱ⁺¹⁾ = x⁽ⁱ⁾ + α⁽ⁱ⁾d⁽ⁱ⁾               ← mise à jour solution
    r⁽ⁱ⁺¹⁾ = r⁽ⁱ⁾ − α⁽ⁱ⁾Ad⁽ⁱ⁾               ← mise à jour résidu
    β⁽ⁱ⁾ = (r⁽ⁱ⁺¹⁾, r⁽ⁱ⁺¹⁾) / (r⁽ⁱ⁾, r⁽ⁱ⁾)  ← coefficient conjugaison
    d⁽ⁱ⁺¹⁾ = r⁽ⁱ⁺¹⁾ + β⁽ⁱ⁾d⁽ⁱ⁾              ← nouvelle direction (A-conjuguée)
Jusqu'à convergence (rayon spectral ρ < 1)
```

Avantage vs. méthodes directes (LU) : complexité **linéaire** contre O(n^(7/3)) pour LU, crucial pour les grands systèmes EF.

La parallélisation distribue les lignes de `A` et les composantes de `x` entre les processus.

---

## Performances mesurées

Expériences réalisées sur le cluster Centaure avec la résolution par gradient conjugué :

**Speed-Up** — évolution quasi-linéaire avec le nombre de cœurs (accélération idéale atteinte) grâce au réseau haute vitesse du Centaure.

**Efficacité** — comprise entre **0.4 et 0.6** : environ 40 à 50 % des ressources sont perdues en communications, synchronisations (`MPI_Barrier`) et déséquilibres de charge.

> **Conclusion** : augmenter le nombre de cœurs ne garantit pas une accélération proportionnelle. Un compromis doit être trouvé entre ressources allouées et gain effectif. L'algorithmique (choix de la méthode de communication) est aussi déterminante que la puissance matérielle.

---

## Auteurs

- **Kevin TONGUE** — École Centrale Lyon ENISE
- **Gaston T. KAMDEM** — École Centrale Lyon ENISE

---

