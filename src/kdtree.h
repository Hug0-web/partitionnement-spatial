#ifndef KDTREE_H
#define KDTREE_H

#include "algos.h"

#define FEUILLE_MAX 8   /* une feuille contient au plus 8 objets */

/* Un nœud couvre les objets idx[debut .. debut+nb). Sa boîte englobe
 * les disques de tous ces objets (pas seulement leurs centres) : c'est
 * ce qui permet de traiter des objets de tailles différentes. */
typedef struct {
    Rect boite;
    int gauche, droite;   /* enfants, -1 pour une feuille */
    int debut, nb;
    int axe;              /* 0 = coupe verticale (x), 1 = horizontale (y) */
    float coupe;          /* coordonnée de la coupe (pour la visualisation) */
} Noeud;

/* Que faire quand les objets ont bougé. */
typedef enum {
    KD_REAJUSTE,      /* recalculer les boîtes, garder l'arbre */
    KD_RECONSTRUIT,   /* refaire l'arbre */
    KD_HYBRIDE        /* réajuster, et refaire l'arbre s'il s'est trop dégradé */
} ModeKd;

typedef struct {
    Noeud *noeuds;
    int nb_noeuds, cap_noeuds;
    int *idx;             /* permutation des indices d'objets */
    int n;
    ModeKd mode;
    float aire_reference; /* somme des aires des feuilles juste après la dernière construction */
    int nb_reconstructions;  /* mode hybride : reconstructions faites depuis la création */
} KdTree;

/* Mode hybride : on reconstruit dès que la somme des aires des feuilles
 * dépasse seuil × sa valeur juste après la dernière construction
 * (1,5 par défaut). Sert à l'étude du seuil. */
void kdtree_fixer_seuil(float seuil);

#endif
