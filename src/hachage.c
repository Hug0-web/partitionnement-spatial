/*
 * Grille hachée (spatial hashing, Teschner et al., 2003).
 *
 * Le monde est découpé en cellules de côté c comme pour la grille, mais
 * sans bornes : la cellule (x, y) est rangée dans le seau
 *     h(x, y) = (x × 73856093  XOR  y × 19349663)  modulo nb_seaux
 * (les deux grands nombres premiers sont ceux de Teschner et al.).
 * Le nombre de seaux est la puissance de 2 au-dessus de 2n : la mémoire
 * est en O(n), quelle que soit la taille du monde. En contrepartie,
 * plusieurs cellules peuvent tomber dans le même seau : chaque test
 * vérifie que l'objet appartient bien à la cellule visitée.
 *
 *   Construction   O(n + e)   (e inscriptions, comme la grille)
 *   Zone, voisins, collisions : comme la grille uniforme
 *   Mise à jour    O(n + e)   : on reconstruit
 *   Mémoire        8 × seaux + 4 e octets, soit 22 à 32 octets par objet
 *                  (selon l'arrondi du nombre de seaux à une puissance de 2)
 *
 * Taille automatique des cellules : 10 fois le rayon moyen des objets.
 * Elle ne dépend pas de la taille du monde, et l'étude de la grille
 * (section 5.5 du rapport) place l'optimum entre 10 et 20 rayons.
 */
#include "hachage.h"

#include <stdlib.h>
#include <string.h>

static float taille_imposee = 0.0f;
static int seaux_imposes = 0;

void hachage_fixer(float taille, int nb_seaux)
{
    taille_imposee = taille;
    seaux_imposes = nb_seaux;
}

/* Cellule (sans bornes) qui contient la coordonnée v. */
static inline int cellule(const Hachage *h, float v)
{
    return (int)floorf(v * h->inv);
}

static inline uint32_t seau(const Hachage *h, int x, int y)
{
    return (((uint32_t)x * 73856093u) ^ ((uint32_t)y * 19349663u)) & h->masque;
}

/* Remplit la table par dénombrement, comme la grille. Un objet qui
 * couvre deux cellules du même seau n'y est inscrit qu'une fois
 * (tableau « dernier »), sinon il serait compté deux fois. */
static void remplir(Hachage *h, const Monde *m)
{
    size_t nb_seaux = (size_t)h->masque + 1;
    int *debut = h->debut;
    memset(debut, 0, (nb_seaux + 1) * sizeof *debut);
    memset(h->dernier, 0, nb_seaux * sizeof *h->dernier);

    /* Passe 1 : compter, et noter l'étendue des cellules occupées. */
    size_t total = 0;
    h->cx_min = h->cy_min = 0;
    h->cx_max = h->cy_max = -1;
    for (int i = 0; i < m->n; i++) {
        Rect b = boite_objet(&m->obj[i]);
        int x0 = cellule(h, b.minx), x1 = cellule(h, b.maxx);
        int y0 = cellule(h, b.miny), y1 = cellule(h, b.maxy);
        if (i == 0 || x0 < h->cx_min) h->cx_min = x0;
        if (i == 0 || y0 < h->cy_min) h->cy_min = y0;
        if (i == 0 || x1 > h->cx_max) h->cx_max = x1;
        if (i == 0 || y1 > h->cy_max) h->cy_max = y1;
        for (int y = y0; y <= y1; y++) {
            for (int x = x0; x <= x1; x++) {
                uint32_t s = seau(h, x, y);
                if (h->dernier[s] == i + 1) continue;   /* déjà inscrit dans ce seau */
                h->dernier[s] = i + 1;
                debut[s + 1]++;
                total++;
            }
        }
    }
    if (total > h->cap_contenu) {
        h->contenu = xrealloc(h->contenu, total * sizeof *h->contenu);
        h->cap_contenu = total;
    }
    for (size_t s = 1; s <= nb_seaux; s++) debut[s] += debut[s - 1];

    /* Passe 2 : placer (debut[s] sert de curseur, puis décalage). */
    memset(h->dernier, 0, nb_seaux * sizeof *h->dernier);
    for (int i = 0; i < m->n; i++) {
        Rect b = boite_objet(&m->obj[i]);
        int x0 = cellule(h, b.minx), x1 = cellule(h, b.maxx);
        int y0 = cellule(h, b.miny), y1 = cellule(h, b.maxy);
        for (int y = y0; y <= y1; y++) {
            for (int x = x0; x <= x1; x++) {
                uint32_t s = seau(h, x, y);
                if (h->dernier[s] == i + 1) continue;
                h->dernier[s] = i + 1;
                h->contenu[debut[s]++] = i;
            }
        }
    }
    memmove(debut + 1, debut, nb_seaux * sizeof *debut);
    debut[0] = 0;
    h->nb_entrees = total;
}

static void *hachage_construire(const Monde *m)
{
    Hachage *h = xcalloc(1, sizeof *h);
    float c = taille_imposee;
    if (c <= 0) {
        double somme = 0;
        for (int i = 0; i < m->n; i++) somme += m->obj[i].r;
        c = 10.0f * (m->n > 0 ? (float)(somme / m->n) : 1.0f);
    }
    h->taille = c;
    h->inv = 1.0f / c;
    size_t nb_seaux = 16;
    size_t voulu = seaux_imposes > 0 ? (size_t)seaux_imposes : 2 * (size_t)m->n;
    while (nb_seaux < voulu) nb_seaux *= 2;
    h->masque = (uint32_t)(nb_seaux - 1);
    h->debut = xmalloc((nb_seaux + 1) * sizeof *h->debut);
    h->dernier = xmalloc(nb_seaux * sizeof *h->dernier);
    remplir(h, m);
    return h;
}

static void hachage_liberer(void *s)
{
    Hachage *h = s;
    free(h->debut);
    free(h->contenu);
    free(h->dernier);
    free(h);
}

/* Parcourt les cellules qui touchent la zone (limitées à celles occupées
 * par des objets) et garde les objets qui la touchent. */
static void parcourir(const Hachage *h, const Monde *m, Rect zone,
                      int cercle, float cx, float cy, float rayon, ListeInt *res)
{
    int x0 = cellule(h, zone.minx), x1 = cellule(h, zone.maxx);
    int y0 = cellule(h, zone.miny), y1 = cellule(h, zone.maxy);
    if (x0 < h->cx_min) x0 = h->cx_min;
    if (y0 < h->cy_min) y0 = h->cy_min;
    if (x1 > h->cx_max) x1 = h->cx_max;
    if (y1 > h->cy_max) y1 = h->cy_max;
    for (int y = y0; y <= y1; y++) {
        for (int x = x0; x <= x1; x++) {
            uint32_t s = seau(h, x, y);
            for (int e = h->debut[s]; e < h->debut[s + 1]; e++) {
                int i = h->contenu[e];
                const Objet *o = &m->obj[i];
                Rect b = boite_objet(o);
                if (!rects_se_coupent(b, zone)) continue;
                /* Comme pour la grille : l'objet n'est compté que dans la
                 * cellule du coin bas-gauche de l'intersection objet ∩ zone.
                 * Ce test écarte aussi les objets d'autres cellules rangées
                 * dans le même seau (ce coin est forcément dans leur boîte). */
                if (cellule(h, fmaxf(b.minx, zone.minx)) != x
                    || cellule(h, fmaxf(b.miny, zone.miny)) != y) continue;
                int dedans = cercle ? objet_dans_cercle(o, cx, cy, rayon)
                                    : objet_dans_rect(o, zone);
                if (dedans) liste_ajouter(res, i);
            }
        }
    }
}

static void hachage_rect(void *s, const Monde *m, Rect zone, ListeInt *res)
{
    parcourir(s, m, zone, 0, 0, 0, 0, res);
}

static void hachage_cercle(void *s, const Monde *m, float cx, float cy, float rayon, ListeInt *res)
{
    parcourir(s, m, boite_cercle(cx, cy, rayon), 1, cx, cy, rayon, res);
}

/* Propose au tas les objets dont le centre est dans la cellule (x, y). */
static void visiter(const Hachage *h, const Monde *m, int x, int y, float px, float py, TasK *tas)
{
    uint32_t s = seau(h, x, y);
    for (int e = h->debut[s]; e < h->debut[s + 1]; e++) {
        int i = h->contenu[e];
        const Objet *o = &m->obj[i];
        if (cellule(h, o->x) != x || cellule(h, o->y) != y) continue;
        tas_proposer(tas, i, dist2(px, py, o->x, o->y));
    }
}

/* k voisins : anneaux de cellules autour du point, comme la grille, en
 * se limitant au rectangle des cellules occupées. */
static int hachage_kppv(void *s, const Monde *m, float px, float py, int k, int *ids, float *d2)
{
    const Hachage *h = s;
    if (k <= 0 || m->n == 0) return 0;
    TasK tas;
    tas_init(&tas, k, ids, d2);
    int qx = cellule(h, px), qy = cellule(h, py);
    if (qx < h->cx_min) qx = h->cx_min;
    if (qx > h->cx_max) qx = h->cx_max;
    if (qy < h->cy_min) qy = h->cy_min;
    if (qy > h->cy_max) qy = h->cy_max;

    for (int r = 0;; r++) {
        for (int y = qy - r; y <= qy + r; y++) {
            if (y < h->cy_min || y > h->cy_max) continue;
            if (y == qy - r || y == qy + r) {
                int x0 = qx - r < h->cx_min ? h->cx_min : qx - r;
                int x1 = qx + r > h->cx_max ? h->cx_max : qx + r;
                for (int x = x0; x <= x1; x++) visiter(h, m, x, y, px, py, &tas);
            } else {
                if (qx - r >= h->cx_min) visiter(h, m, qx - r, y, px, py, &tas);
                if (qx + r <= h->cx_max) visiter(h, m, qx + r, y, px, py, &tas);
            }
        }
        float borne = INFINITY;
        int reste = 0;
        if (qx - r > h->cx_min) { reste = 1; borne = fminf(borne, px - (float)(qx - r) * h->taille); }
        if (qx + r < h->cx_max) { reste = 1; borne = fminf(borne, (float)(qx + r + 1) * h->taille - px); }
        if (qy - r > h->cy_min) { reste = 1; borne = fminf(borne, py - (float)(qy - r) * h->taille); }
        if (qy + r < h->cy_max) { reste = 1; borne = fminf(borne, (float)(qy + r + 1) * h->taille - py); }
        if (!reste) break;
        float lim = borne - 0.001f * h->taille;
        if (tas.n == k && lim > 0 && tas_seuil(&tas) <= lim * lim) break;
    }
    return tas_finir(&tas);
}

/* Collisions : on parcourt les seaux. Une paire n'est comptée que dans
 * le seau de la cellule X qui contient le coin bas-gauche de
 * l'intersection des deux boîtes ; les deux objets couvrent X, ils sont
 * donc bien inscrits (une seule fois) dans ce seau.
 * Découpage pour le parallélisme : la part p traite les seaux p, p + nb... */
static void hachage_paires(void *s, const Monde *m, int part, int nb_parts, ListePaires *res)
{
    const Hachage *h = s;
    for (uint32_t sc = (uint32_t)part; sc <= h->masque; sc += (uint32_t)nb_parts) {
        int d = h->debut[sc], f = h->debut[sc + 1];
        for (int ea = d; ea < f; ea++) {
            int i = h->contenu[ea];
            Rect bi = boite_objet(&m->obj[i]);
            for (int eb = ea + 1; eb < f; eb++) {
                int j = h->contenu[eb];
                Rect bj = boite_objet(&m->obj[j]);
                if (!rects_se_coupent(bi, bj)) continue;
                int x = cellule(h, fmaxf(bi.minx, bj.minx));
                int y = cellule(h, fmaxf(bi.miny, bj.miny));
                if (seau(h, x, y) != sc) continue;
                if (collision(&m->obj[i], &m->obj[j])) paires_ajouter(res, i, j);
            }
        }
    }
}

static void hachage_mettre_a_jour(void *s, const Monde *m)
{
    remplir(s, m);
}

static size_t hachage_memoire(const void *s)
{
    const Hachage *h = s;
    size_t nb_seaux = (size_t)h->masque + 1;
    return sizeof *h + (nb_seaux + 1) * sizeof *h->debut + nb_seaux * sizeof *h->dernier
         + h->cap_contenu * sizeof *h->contenu;
}

const Algo ALGO_HACHAGE = {
    "grille-hachee", hachage_construire, hachage_liberer, hachage_rect, hachage_cercle,
    hachage_kppv, hachage_paires, hachage_mettre_a_jour, hachage_memoire
};
