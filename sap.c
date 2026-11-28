/*
 * Sweep and Prune sur un axe (Cohen et al., 1995, I-COLLIDE).
 *
 * On garde les objets triés selon le bord gauche (minx) de leur boîte.
 * Pour trouver les collisions, on balaie la liste : un objet ne peut
 * toucher que les objets suivants dont le bord gauche est avant son bord
 * droit. Quand les objets bougent un peu, la liste reste presque triée :
 * un tri par insertion la remet en ordre en O(n + échanges). C'est la
 * cohérence temporelle exploitée par I-COLLIDE.
 *
 *   Construction   O(n log n)  : tri initial
 *   Zone           O(log n + objets dont l'intervalle en x touche la zone)
 *   k voisins      balayage à gauche et à droite du point, avec arrêt
 *   Collisions     O(n + paires qui se chevauchent en x)
 *   Mise à jour    O(n + échanges), O(n) si les objets bougent peu
 *   Mémoire        8 n octets
 *
 * Point faible : un seul axe. Si beaucoup d'objets ont la même abscisse
 * (amas, gros objets), les intervalles en x se chevauchent tous et on
 * retombe vers O(n²).
 */
#include "algos.h"

#include <stdlib.h>

typedef struct {
    float min;   /* bord gauche de la boîte de l'objet */
    int id;
} Entree;

typedef struct {
    Entree *e;
    int n;
    float largeur_max;   /* plus grande largeur de boîte (2 × rayon max) */
} Sap;

static int comparer_entrees(const void *pa, const void *pb)
{
    const Entree *a = pa, *b = pb;
    if (a->min != b->min) return a->min < b->min ? -1 : 1;
    return (a->id > b->id) - (a->id < b->id);
}

static void calculer_largeur_max(Sap *s, const Monde *m)
{
    float rmax = 0;
    for (int i = 0; i < m->n; i++) rmax = fmaxf(rmax, m->obj[i].r);
    /* Marge pour les arrondis : (x + r) - (x - r) peut dépasser 2r d'un ulp. */
    s->largeur_max = 2 * rmax * 1.0001f + 1e-3f;
}

static void *sap_construire(const Monde *m)
{
    Sap *s = xcalloc(1, sizeof *s);
    s->n = m->n;
    s->e = xmalloc((size_t)(m->n > 0 ? m->n : 1) * sizeof *s->e);
    for (int i = 0; i < m->n; i++) {
        s->e[i].min = boite_objet(&m->obj[i]).minx;
        s->e[i].id = i;
    }
    qsort(s->e, (size_t)s->n, sizeof *s->e, comparer_entrees);
    calculer_largeur_max(s, m);
    return s;
}

static void sap_liberer(void *p)
{
    Sap *s = p;
    free(s->e);
    free(s);
}

/* Première position dont le bord gauche est >= v (recherche dichotomique). */
static int premier_superieur(const Sap *s, float v)
{
    int g = 0, d = s->n;
    while (g < d) {
        int mil = g + (d - g) / 2;
        if (s->e[mil].min < v) g = mil + 1;
        else d = mil;
    }
    return g;
}

/* Un objet qui touche la zone a son bord gauche dans
 * [zone.minx - largeur_max, zone.maxx] : on ne parcourt que cette tranche. */
static void parcourir(const Sap *s, const Monde *m, Rect zone,
                      int cercle, float cx, float cy, float rayon, ListeInt *res)
{
    for (int p = premier_superieur(s, zone.minx - s->largeur_max);
         p < s->n && s->e[p].min <= zone.maxx; p++) {
        int i = s->e[p].id;
        int dedans = cercle ? objet_dans_cercle(&m->obj[i], cx, cy, rayon)
                            : objet_dans_rect(&m->obj[i], zone);
        if (dedans) liste_ajouter(res, i);
    }
}

static void sap_rect(void *s, const Monde *m, Rect zone, ListeInt *res)
{
    parcourir(s, m, zone, 0, 0, 0, 0, res);
}

static void sap_cercle(void *s, const Monde *m, float cx, float cy, float rayon, ListeInt *res)
{
    parcourir(s, m, boite_cercle(cx, cy, rayon), 1, cx, cy, rayon, res);
}

/* k voisins : on part de la position du point dans la liste et on
 * s'éloigne des deux côtés. Le centre d'un objet est entre son bord
 * gauche et son bord gauche + largeur_max : dès que cet intervalle est
 * trop loin en x, on arrête ce côté. */
static int sap_kppv(void *p, const Monde *m, float x, float y, int k, int *ids, float *d2)
{
    const Sap *s = p;
    if (k <= 0 || m->n == 0) return 0;
    TasK tas;
    tas_init(&tas, k, ids, d2);
    int depart = premier_superieur(s, x);
    for (int q = depart; q < s->n; q++) {             /* vers la droite */
        float ecart = s->e[q].min - x;                  /* >= 0 */
        if (ecart * ecart >= tas_seuil(&tas)) break;
        const Objet *o = &m->obj[s->e[q].id];
        tas_proposer(&tas, s->e[q].id, dist2(x, y, o->x, o->y));
    }
    for (int q = depart - 1; q >= 0; q--) {           /* vers la gauche */
        float ecart = x - (s->e[q].min + s->largeur_max);
        if (ecart > 0 && ecart * ecart >= tas_seuil(&tas)) break;
        const Objet *o = &m->obj[s->e[q].id];
        tas_proposer(&tas, s->e[q].id, dist2(x, y, o->x, o->y));
    }
    return tas_finir(&tas);
}

/* Découpage pour le parallélisme : la part p traite les positions
 * a = p, p + nb, p + 2 nb... de la liste triée. */
static void sap_paires(void *p, const Monde *m, int part, int nb_parts, ListePaires *res)
{
    const Sap *s = p;
    for (int a = part; a < s->n; a += nb_parts) {
        int i = s->e[a].id;
        const Objet *oi = &m->obj[i];
        float maxx = boite_objet(oi).maxx;
        /* Les suivants commencent après oi ; dès que l'un commence après
         * le bord droit de oi, tous les autres aussi. */
        for (int b = a + 1; b < s->n && s->e[b].min <= maxx; b++) {
            int j = s->e[b].id;
            if (collision(oi, &m->obj[j])) paires_ajouter(res, i, j);
        }
    }
}

static void sap_mettre_a_jour(void *p, const Monde *m)
{
    Sap *s = p;
    for (int q = 0; q < s->n; q++) s->e[q].min = boite_objet(&m->obj[s->e[q].id]).minx;
    /* Tri par insertion : O(n) si chaque objet n'a doublé que quelques voisins. */
    for (int q = 1; q < s->n; q++) {
        Entree v = s->e[q];
        int r = q - 1;
        while (r >= 0 && comparer_entrees(&s->e[r], &v) > 0) {
            s->e[r + 1] = s->e[r];
            r--;
        }
        s->e[r + 1] = v;
    }
    calculer_largeur_max(s, m);
}

static size_t sap_memoire(const void *p)
{
    const Sap *s = p;
    return sizeof *s + (size_t)(s->n > 0 ? s->n : 1) * sizeof *s->e;
}

const Algo ALGO_SAP = {
    "sweep-and-prune", sap_construire, sap_liberer, sap_rect, sap_cercle,
    sap_kppv, sap_paires, sap_mettre_a_jour, sap_memoire
};
