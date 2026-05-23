# ⚡ Calculs Intensifs — MPI en C++

Travaux pratiques de calcul parallèle réalisés sur le cluster **Centaure** de l'École des Mines de Saint-Étienne (656 cœurs). On a codé en C++ avec MPI, progressé exercice par exercice, et fini par résoudre un vrai système linéaire éléments finis en parallèle.

> Projet réalisé en binôme avec **Gaston T. KAMDEM** — École Centrale Lyon ENISE

---

## Ce qu'on a fait

### 👋 Hello World — premiers pas avec MPI

Avant de faire quoi que ce soit, il faut comprendre comment MPI fonctionne. Chaque processus a un **rang** (son identifiant), et ils communiquent tous dans un espace partagé appelé `MPI_COMM_WORLD`. C'est le b.a.-ba pour tout ce qui suit.

---

### ➕ Somme des n premiers entiers — le vrai départ

L'idée : calculer `1 + 2 + ... + n` mais en le répartissant sur plusieurs cœurs. Chaque cœur fait sa portion, puis on agrège tout.

Ce qui est intéressant ici, c'est qu'on a comparé **3 façons de faire communiquer les processus** entre eux :

| Méthode | Comment ça marche | Coût |
|---|---|---|
| 🐌 Naïve | Tout le monde parle à tout le monde | O(p²) |
| 🚶 Intermédiaire | Tout le monde parle au processus 0 | O(p) |
| 🚀 Optimale | Réduction par arbre binaire (pair/impair) | O(log p) |

Morale : **l'algo compte autant que le matériel**.

---

### 🔢 Gradient Conjugué — le gros morceau

On résout un système `Ax = b` comme ceux qu'on rencontre en éléments finis. La matrice `A` est grande, creuse, symétrique définie positive — exactement le genre de problème où les méthodes directes (type LU) sont trop lourdes.

Le gradient conjugué, c'est une méthode **itérative** : on part d'une solution approchée et on converge vers la bonne en suivant des directions conjuguées. Complexité linéaire, bien adapté au calcul distribué.

On a analysé les performances selon le nombre de cœurs : le **speed-up est quasi-linéaire** (accélération idéale ✅), mais l'**efficacité plafonne entre 0.4 et 0.6** — 40 à 50% des ressources partent en fumée dans les communications et synchronisations. Pas si surprenant, c'est le prix du parallélisme.

---

## Structure du projet

```
CalculsIntensifs/
├── hello_exercice/     # Hello World MPI
├── somme_exercice/     # Somme parallèle + 3 stratégies de comm
└── GC_complet/         # Gradient Conjugué parallèle
```

---

## Lancer le code

```bash
# Compiler
cd somme_exercice && make

# Exécuter sur 4 cœurs
mpirun -np 4 ./somme

# Sur le cluster
mpirun -np 16 ./GC_complet/gradient_conjugue
```

---

## Ce qu'on a retenu

Ajouter des cœurs ne suffit pas. Ce qui fait la différence, c'est la façon dont les processus communiquent — trop de messages, et le gain s'évapore. Sur ce projet, l'algorithmique a été au moins aussi importante que la puissance brute du cluster.

---

*C++ · MPI · HPC · École Centrale Lyon ENISE*
