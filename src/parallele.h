#ifndef PARALLELE_H
#define PARALLELE_H

#include "algos.h"

#define FILS_MAX 64

/* Toutes les paires en collision, calcul réparti sur nb_fils fils
 * d'exécution (threads POSIX). Même résultat que a->paires(s, m, 0, 1, res),
 * dans un autre ordre. */
void paires_paralleles(const Algo *a, void *s, const Monde *m, int nb_fils, ListePaires *res);

/* Chaque objet cherche ses k plus proches voisins (comme chaque ennemi
 * d'un jeu qui regarde autour de lui). Les voisins de l'objet i sont
 * écrits dans ids[i*k .. i*k + k). Calcul réparti sur nb_fils fils. */
void kppv_tous_paralleles(const Algo *a, void *s, const Monde *m, int k, int nb_fils, int *ids);

#endif
