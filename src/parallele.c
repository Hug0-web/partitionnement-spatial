/*
 * Version parallélisée (threads POSIX).
 *
 * Les structures ne sont que lues pendant les requêtes et la recherche
 * des collisions (aucune variable modifiée, pas de tableau de marquage) :
 * plusieurs fils peuvent donc les interroger en même temps sans verrou.
 * Chaque fil écrit dans ses propres résultats, fusionnés à la fin.
 *
 * Collisions : chaque structure sait calculer la part p d'un découpage
 * en nb morceaux disjoints (paramètres part / nb_parts de a->paires) ;
 * le fil p calcule la part p.
 * Voisins de tous les objets : le fil p traite les objets p, p + nb, ...
 *
 * Pourquoi un découpage en tourniquet plutôt qu'en blocs contigus :
 * dans une scène en amas, les objets proches sont aussi voisins dans les
 * structures (lignes de la grille, seaux, feuilles du kd-tree, positions
 * du balayage). Avec des blocs, un fil pourrait recevoir tout un amas
 * dense et travailler seul pendant que les autres attendent. En prenant
 * un élément sur nb, chaque fil reçoit un peu de chaque zone, et la
 * charge reste à peu près égale sans avoir à la mesurer.
 */
#define _POSIX_C_SOURCE 200112L   /* pthread */

#include "parallele.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const Algo *a;
    void *s;
    const Monde *m;
    int part, nb_parts;
    int k;              /* voisins : nombre demandé */
    int *ids;           /* voisins : tableau de sortie commun (zones disjointes) */
    ListePaires paires; /* collisions : résultats propres à ce fil */
} Tache;

static void *travail_paires(void *arg)
{
    Tache *t = arg;
    t->a->paires(t->s, t->m, t->part, t->nb_parts, &t->paires);
    return NULL;
}

static void *travail_kppv(void *arg)
{
    Tache *t = arg;
    int *voisins = xmalloc((size_t)t->k * sizeof *voisins);
    float *d2 = xmalloc((size_t)t->k * sizeof *d2);
    /* L'objet i écrit dans ids[i*k .. i*k + k - 1] : ces zones ne se
     * chevauchent pas d'un fil à l'autre, d'où l'absence de verrou. */
    for (int i = t->part; i < t->m->n; i += t->nb_parts) {
        const Objet *o = &t->m->obj[i];
        int nb = t->a->kppv(t->s, t->m, o->x, o->y, t->k, voisins, d2);
        int *sortie = t->ids + (size_t)i * t->k;
        for (int v = 0; v < t->k; v++) sortie[v] = v < nb ? voisins[v] : -1;
    }
    free(voisins);
    free(d2);
    return NULL;
}

/* Lance fonction sur les nb tâches : le fil appelant fait la tâche 0,
 * des fils créés font les autres, puis on attend qu'ils finissent. */
static void executer(void *(*fonction)(void *), Tache *taches, int nb)
{
    pthread_t fils[FILS_MAX];
    int lance[FILS_MAX] = {0};
    for (int p = 1; p < nb; p++) {
        /* Si le fil ne peut pas être créé, on fait le travail ici. */
        if (pthread_create(&fils[p], NULL, fonction, &taches[p]) == 0) lance[p] = 1;
        else fonction(&taches[p]);
    }
    fonction(&taches[0]);
    for (int p = 1; p < nb; p++) {
        if (lance[p]) pthread_join(fils[p], NULL);
    }
}

static int borner_fils(int nb_fils)
{
    if (nb_fils < 1) return 1;
    return nb_fils > FILS_MAX ? FILS_MAX : nb_fils;
}

void paires_paralleles(const Algo *a, void *s, const Monde *m, int nb_fils, ListePaires *res)
{
    nb_fils = borner_fils(nb_fils);
    Tache taches[FILS_MAX];
    for (int p = 0; p < nb_fils; p++) {
        taches[p] = (Tache){ a, s, m, p, nb_fils, 0, NULL, {0} };
    }
    executer(travail_paires, taches, nb_fils);
    /* Fusion des résultats de chaque fil, dans l'ordre des parts : on
     * copie après coup plutôt que d'ajouter à une liste commune pendant
     * le calcul, qui demanderait un verrou à chaque paire trouvée. */
    for (int p = 0; p < nb_fils; p++) {
        ListePaires *l = &taches[p].paires;
        if (res->n + l->n > res->cap) {
            res->cap = res->n + l->n;
            res->v = xrealloc(res->v, res->cap * sizeof *res->v);
        }
        if (l->n > 0) memcpy(res->v + res->n, l->v, l->n * sizeof *l->v);
        res->n += l->n;
        paires_liberer(l);
    }
}

void kppv_tous_paralleles(const Algo *a, void *s, const Monde *m, int k, int nb_fils, int *ids)
{
    if (k <= 0) return;
    nb_fils = borner_fils(nb_fils);
    Tache taches[FILS_MAX];
    for (int p = 0; p < nb_fils; p++) {
        taches[p] = (Tache){ a, s, m, p, nb_fils, k, ids, {0} };
    }
    executer(travail_kppv, taches, nb_fils);
}
