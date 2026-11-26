#ifndef GRILLE_H
#define GRILLE_H

#include "algos.h"

/* Grille uniforme stockée en « lignes compressées » : les objets de la
 * cellule c sont contenu[debut[c] .. debut[c+1]). Un objet est inscrit
 * dans toutes les cellules que touche sa boîte englobante. */
typedef struct {
    float taille, inv;      /* côté d'une cellule et son inverse */
    int nx, ny;             /* nombre de cellules en x et en y */
    int *debut;             /* nx*ny + 1 entrées */
    int *contenu;
    size_t nb_entrees;      /* nombre d'inscriptions (>= n) */
    size_t cap_contenu;
} Grille;

/* Côté des cellules imposé pour les prochaines constructions
 * (0 = choix automatique). Sert à l'étude de sensibilité. */
void grille_fixer_taille(float taille);

#endif
