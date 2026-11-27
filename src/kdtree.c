/*
 * kd-tree (Bentley, 1975) sur les centres des objets.
 *
 * Chaque nœud coupe ses objets en deux moitiés égales selon la médiane
 * de l'axe où les centres sont le plus étalés. On s'arrête à 8 objets
 * par feuille. Chaque nœud garde la boîte qui englobe les disques de
 * ses objets : un objet n'est rangé qu'une fois, même s'il est gros, et
 * les requêtes élaguent les sous-arbres dont la boîte ne touche pas la
 * zone. (Garder ces boîtes rapproche ce kd-tree d'un BVH.)
 *
 *   Construction   O(n log n)  : sélection de la médiane en O(n) par niveau
 *   Zone           O(√n + m) en moyenne pour une petite zone (m résultats)
 *   k voisins      O(log n + k) en moyenne
 *   Collisions     O(n log n + paires) : une requête par feuille
 *   Mise à jour    réajustement O(n), ou reconstruction O(n log n)
 *   Mémoire        n/2 nœuds de 40 octets + 4 n octets
 *
 * Trois variantes de mise à jour quand les objets bougent :
 *   - réajustée    : on garde l'arbre et on recalcule seulement les
 *                    boîtes, des feuilles vers la racine, en O(n). Rapide,
 *                    mais les boîtes grossissent et se chevauchent au fil
 *                    des déplacements : les requêtes ralentissent ;
 *   - reconstruite : on refait l'arbre en O(n log n) à chaque pas ;
 *   - hybride      : on réajuste, et on reconstruit quand la somme des
 *                    aires des feuilles a trop grossi. Le coût d'une
 *                    requête est à peu près proportionnel à cette somme
 *                    (c'est l'idée de l'heuristique des surfaces, SAH,
 *                    MacDonald et Booth, 1990) : elle mesure la dégradation.
 */
#include "kdtree.h"

#include <stdlib.h>

#define PILE_MAX 128   /* profondeur max ≈ log2(n/4) + 1, largement suffisant */

static float seuil_hybride = 1.5f;

void kdtree_fixer_seuil(float seuil)
{
    seuil_hybride = seuil;
}

static inline float coord(const Objet *o, int axe)
{
    return axe == 0 ? o->x : o->y;
}

/* Sélection rapide (Hoare) : réordonne idx[0 .. n) pour que idx[k] soit
 * l'élément de rang k selon l'axe, les plus petits avant, les plus
 * grands après. O(n) en moyenne. */
static void selection(int *idx, int n, int k, const Objet *obj, int axe)
{
    int g = 0, d = n - 1;
    while (g < d) {
        float pivot = coord(&obj[idx[(g + d) / 2]], axe);
        int i = g, j = d;
        while (i <= j) {
            while (coord(&obj[idx[i]], axe) < pivot) i++;
            while (coord(&obj[idx[j]], axe) > pivot) j--;
            if (i <= j) {
                int t = idx[i]; idx[i] = idx[j]; idx[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j) d = j;          /* le rang k est à gauche */
        else if (k >= i) g = i;     /* à droite */
        else break;                 /* entre les deux : égal au pivot */
    }
}

static Rect boite_feuille(const KdTree *t, const Noeud *nd, const Objet *obj)
{
    Rect b = rect_vide();
    for (int e = nd->debut; e < nd->debut + nd->nb; e++) {
        b = rect_union(b, boite_objet(&obj[t->idx[e]]));
    }
    return b;
}

/* Construit le sous-arbre des objets idx[debut .. debut+nb) et renvoie
 * l'indice de sa racine. Les nœuds sont numérotés dans l'ordre préfixe :
 * un enfant a toujours un indice plus grand que son parent. */
static int construire_rec(KdTree *t, const Objet *obj, int debut, int nb)
{
    int id = t->nb_noeuds++;
    Noeud *nd = &t->noeuds[id];
    nd->debut = debut;
    nd->nb = nb;
    nd->gauche = nd->droite = -1;
    nd->axe = 0;
    nd->coupe = 0;

    if (nb <= FEUILLE_MAX) {
        nd->boite = boite_feuille(t, nd, obj);
        return id;
    }

    /* Axe de coupe : celui où les centres sont le plus étalés. */
    Rect c = rect_vide();
    for (int e = debut; e < debut + nb; e++) {
        const Objet *o = &obj[t->idx[e]];
        Rect p = { o->x, o->y, o->x, o->y };
        c = rect_union(c, p);
    }
    int axe = (c.maxx - c.minx) >= (c.maxy - c.miny) ? 0 : 1;
    int mil = nb / 2;
    selection(t->idx + debut, nb, mil, obj, axe);

    int g = construire_rec(t, obj, debut, mil);
    int d = construire_rec(t, obj, debut + mil, nb - mil);
    nd = &t->noeuds[id];
    nd->gauche = g;
    nd->droite = d;
    nd->axe = axe;
    nd->coupe = coord(&obj[t->idx[debut + mil]], axe);
    nd->boite = rect_union(t->noeuds[g].boite, t->noeuds[d].boite);
    return id;
}

/* Somme des aires des boîtes des feuilles : mesure de la qualité de
 * l'arbre (plus elle est grande, plus les requêtes visitent de feuilles). */
static float somme_aires_feuilles(const KdTree *t)
{
    double somme = 0;
    for (int id = 0; id < t->nb_noeuds; id++) {
        const Noeud *nd = &t->noeuds[id];
        if (nd->gauche < 0 && nd->nb > 0)
            somme += (double)(nd->boite.maxx - nd->boite.minx) * (nd->boite.maxy - nd->boite.miny);
    }
    return (float)somme;
}

static void reconstruire(KdTree *t, const Monde *m)
{
    t->nb_noeuds = 0;
    construire_rec(t, m->obj, 0, m->n);
    t->aire_reference = somme_aires_feuilles(t);
}

static void *kd_construire_mode(const Monde *m, ModeKd mode)
{
    KdTree *t = xcalloc(1, sizeof *t);
    t->n = m->n;
    t->mode = mode;
    t->idx = xmalloc((size_t)(m->n > 0 ? m->n : 1) * sizeof *t->idx);
    for (int i = 0; i < m->n; i++) t->idx[i] = i;
    /* Toute feuille issue d'une coupe contient au moins 4 objets (moitié
     * de 9 ou plus), donc au plus n/4 feuilles et n/2 nœuds. */
    t->cap_noeuds = m->n / 2 + 2;
    t->noeuds = xmalloc((size_t)t->cap_noeuds * sizeof *t->noeuds);
    reconstruire(t, m);
    return t;
}

static void *kd_construire_reajuste(const Monde *m)    { return kd_construire_mode(m, KD_REAJUSTE); }
static void *kd_construire_reconstruit(const Monde *m) { return kd_construire_mode(m, KD_RECONSTRUIT); }
static void *kd_construire_hybride(const Monde *m)     { return kd_construire_mode(m, KD_HYBRIDE); }

static void kd_liberer(void *s)
{
    KdTree *t = s;
    free(t->noeuds);
    free(t->idx);
    free(t);
}

/* Parcours en profondeur avec une pile : on ne descend que dans les
 * nœuds dont la boîte touche la zone. */
static void parcourir(const KdTree *t, const Monde *m, Rect zone,
                      int cercle, float cx, float cy, float rayon, ListeInt *res)
{
    if (t->nb_noeuds == 0) return;
    int pile[PILE_MAX], sp = 0;
    pile[sp++] = 0;
    while (sp > 0) {
        const Noeud *nd = &t->noeuds[pile[--sp]];
        if (!rects_se_coupent(nd->boite, zone)) continue;
        if (nd->gauche < 0) {
            for (int e = nd->debut; e < nd->debut + nd->nb; e++) {
                int i = t->idx[e];
                int dedans = cercle ? objet_dans_cercle(&m->obj[i], cx, cy, rayon)
                                    : objet_dans_rect(&m->obj[i], zone);
                if (dedans) liste_ajouter(res, i);
            }
        } else {
            pile[sp++] = nd->gauche;
            pile[sp++] = nd->droite;
        }
    }
}

static void kd_rect(void *s, const Monde *m, Rect zone, ListeInt *res)
{
    parcourir(s, m, zone, 0, 0, 0, 0, res);
}

static void kd_cercle(void *s, const Monde *m, float cx, float cy, float rayon, ListeInt *res)
{
    parcourir(s, m, boite_cercle(cx, cy, rayon), 1, cx, cy, rayon, res);
}

/* Descente récursive, enfant le plus proche d'abord ; on abandonne un
 * sous-arbre dès que sa boîte est plus loin que le k-ième candidat. La
 * boîte contient tous les centres : sa distance est un minorant. */
static void kppv_rec(const KdTree *t, const Monde *m, int ni, float x, float y, TasK *tas)
{
    const Noeud *nd = &t->noeuds[ni];
    if (dist2_point_rect(x, y, nd->boite) >= tas_seuil(tas)) return;
    if (nd->gauche < 0) {
        for (int e = nd->debut; e < nd->debut + nd->nb; e++) {
            int i = t->idx[e];
            tas_proposer(tas, i, dist2(x, y, m->obj[i].x, m->obj[i].y));
        }
        return;
    }
    float dg = dist2_point_rect(x, y, t->noeuds[nd->gauche].boite);
    float dd = dist2_point_rect(x, y, t->noeuds[nd->droite].boite);
    int premier = dg <= dd ? nd->gauche : nd->droite;
    int second = dg <= dd ? nd->droite : nd->gauche;
    kppv_rec(t, m, premier, x, y, tas);
    kppv_rec(t, m, second, x, y, tas);
}

static int kd_kppv(void *s, const Monde *m, float x, float y, int k, int *ids, float *d2)
{
    const KdTree *t = s;
    if (k <= 0 || m->n == 0) return 0;
    TasK tas;
    tas_init(&tas, k, ids, d2);
    kppv_rec(t, m, 0, x, y, &tas);
    return tas_finir(&tas);
}

/* Teste toutes les paires entre la feuille a et la feuille b
 * (ou à l'intérieur de a si a == b). */
static void paires_feuilles(const KdTree *t, const Monde *m, const Noeud *a, const Noeud *b,
                            ListePaires *res)
{
    for (int ea = a->debut; ea < a->debut + a->nb; ea++) {
        int i = t->idx[ea];
        const Objet *oi = &m->obj[i];
        /* Dans la même feuille, on ne prend que les objets placés après. */
        int eb0 = a == b ? ea + 1 : b->debut;
        for (int eb = eb0; eb < b->debut + b->nb; eb++) {
            int j = t->idx[eb];
            if (collision(oi, &m->obj[j])) paires_ajouter(res, i, j);
        }
    }
}

/* Collisions : chaque feuille interroge l'arbre avec sa boîte et teste
 * ses objets contre ceux des feuilles trouvées. Deux objets en collision
 * ont des boîtes qui se touchent, donc leurs feuilles aussi : aucune
 * paire n'est oubliée. Une paire de feuilles n'est traitée que depuis
 * celle de plus petit numéro, pour ne compter chaque paire qu'une fois.
 * C'est 8 fois moins de descentes dans l'arbre qu'une requête par objet.
 * Mesuré à 50 000 objets uniformes : 2,3 ms, contre 6,8 ms avec une
 * requête par objet (et 10,9 ms en parcourant les objets dans leur ordre
 * de création plutôt que dans l'ordre des feuilles, à cause du cache).
 * Découpage pour le parallélisme : la part p traite les feuilles numéro
 * p, p + nb, p + 2 nb... */
static void kd_paires(void *s, const Monde *m, int part, int nb_parts, ListePaires *res)
{
    const KdTree *t = s;
    int pile[PILE_MAX];
    int num_feuille = 0;
    for (int fa = 0; fa < t->nb_noeuds; fa++) {
        const Noeud *a = &t->noeuds[fa];
        if (a->gauche >= 0) continue;            /* seulement les feuilles */
        if (num_feuille++ % nb_parts != part) continue;
        int sp = 0;
        pile[sp++] = 0;
        while (sp > 0) {
            int ni = pile[--sp];
            const Noeud *nd = &t->noeuds[ni];
            if (!rects_se_coupent(nd->boite, a->boite)) continue;
            if (nd->gauche < 0) {
                if (ni >= fa) paires_feuilles(t, m, a, nd, res);
            } else {
                pile[sp++] = nd->gauche;
                pile[sp++] = nd->droite;
            }
        }
    }
}

static void kd_mettre_a_jour(void *s, const Monde *m)
{
    KdTree *t = s;
    if (t->mode == KD_RECONSTRUIT) {
        reconstruire(t, m);
        return;
    }
    /* Réajustement : les enfants ont des indices plus grands que leur
     * parent, un parcours à l'envers traite donc les enfants d'abord. */
    for (int id = t->nb_noeuds - 1; id >= 0; id--) {
        Noeud *nd = &t->noeuds[id];
        if (nd->gauche < 0) nd->boite = boite_feuille(t, nd, m->obj);
        else nd->boite = rect_union(t->noeuds[nd->gauche].boite, t->noeuds[nd->droite].boite);
    }
    /* Hybride : si les feuilles ont trop grossi, on refait l'arbre. */
    if (t->mode == KD_HYBRIDE && somme_aires_feuilles(t) > seuil_hybride * t->aire_reference) {
        reconstruire(t, m);
        t->nb_reconstructions++;
    }
}

static size_t kd_memoire(const void *s)
{
    const KdTree *t = s;
    return sizeof *t + (size_t)t->cap_noeuds * sizeof *t->noeuds
         + (size_t)(t->n > 0 ? t->n : 1) * sizeof *t->idx;
}

const Algo ALGO_KDTREE_REAJUSTE = {
    "kdtree-reajuste", kd_construire_reajuste, kd_liberer, kd_rect, kd_cercle,
    kd_kppv, kd_paires, kd_mettre_a_jour, kd_memoire
};

const Algo ALGO_KDTREE_RECONSTRUIT = {
    "kdtree-reconstruit", kd_construire_reconstruit, kd_liberer, kd_rect, kd_cercle,
    kd_kppv, kd_paires, kd_mettre_a_jour, kd_memoire
};

const Algo ALGO_KDTREE_HYBRIDE = {
    "kdtree-hybride", kd_construire_hybride, kd_liberer, kd_rect, kd_cercle,
    kd_kppv, kd_paires, kd_mettre_a_jour, kd_memoire
};
