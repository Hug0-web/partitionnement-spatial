#ifndef HACHAGE_H
#define HACHAGE_H

#include "algos.h"

/* Grille hachée : même découpage en cellules que la grille uniforme,
 * mais les cellules sont rangées dans une table de « seaux » de taille
 * fixe (puissance de 2, au moins 2n) par une fonction de hachage. La
 * mémoire ne dépend plus de la surface du monde, seulement de n. */
typedef struct {
    float taille, inv;        /* côté d'une cellule et son inverse */
    uint32_t masque;          /* nombre de seaux - 1 */
    int *debut;               /* seaux + 1 entrées (comme la grille) */
    int *contenu;
    int *dernier;             /* par seau : dernier objet inscrit + 1 (anti-doublon) */
    size_t nb_entrees, cap_contenu;
    int cx_min, cy_min, cx_max, cy_max;  /* cellules occupées par les objets */
} Hachage;

/* Côté des cellules et nombre de seaux imposés pour les prochaines
 * constructions (0 = choix automatique). Sert aux tests et aux mesures. */
void hachage_fixer(float taille, int nb_seaux);

#endif
