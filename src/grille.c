/*
 * Grille uniforme (Teschner et al., 2003, pour la variante hachée).
 *
 * Le monde est découpé en cellules carrées de côté c. Chaque objet est
 * inscrit dans les cellules que touche sa boîte englobante. Une requête
 * ne regarde que les cellules qui touchent la zone.
 *
 * Avec C cellules, n objets et e inscriptions (e = n si les objets sont
 * petits devant c) :
 *   Construction   O(n + C + e)  : tri par dénombrement, deux passes
 *   Zone           O(cellules touchées + objets dans ces cellules)
 *   k voisins      O(anneaux de cellules parcourus), O(k) en moyenne
 *   Collisions     O(n + somme des m² par cellule), O(n) si la densité est bornée
 *   Mise à jour    O(n + C + e)  : on reconstruit (c'est le plus rapide ici)
 *   Mémoire        4 (C + 1) + 4 e octets
 *
 * Choix automatique de c : la distance moyenne entre objets
 * (environ un objet par cellule), au moins deux fois le rayon moyen.
 */
#include "grille.h"

#include <stdlib.h>
#include <string.h>

static float taille_imposee = 0.0f;

void grille_fixer_taille(float taille)
{
    taille_imposee = taille;
}

/* Cellule contenant la coordonnée x, ramenée dans la grille.
 * La fonction est croissante : c'est ce qui garantit la cohérence entre
 * l'inscription des objets et leur recherche. */
static inline int cellule_x(const Grille *g, float x)
{
    float f = x * g->inv;
    if (f < 0) return 0;
    if (f >= (float)g->nx) return g->nx - 1;
    return (int)f;
}

static inline int cellule_y(const Grille *g, float y)
{
    float f = y * g->inv;
    if (f < 0) return 0;
    if (f >= (float)g->ny) return g->ny - 1;
    return (int)f;
}

/* Remplit la grille par dénombrement :
 * 1. compter les inscriptions de chaque cellule ;
 * 2. sommes cumulées -> début de chaque cellule ;
 * 3. placer les objets. */
static void remplir(Grille *g, const Monde *m)
{
    size_t nc = (size_t)g->nx * (size_t)g->ny;
    int *debut = g->debut;
    memset(debut, 0, (nc + 1) * sizeof *debut);

    size_t total = 0;
    for (int i = 0; i < m->n; i++) {
        Rect b = boite_objet(&m->obj[i]);
        int x0 = cellule_x(g, b.minx), x1 = cellule_x(g, b.maxx);
        int y0 = cellule_y(g, b.miny), y1 = cellule_y(g, b.maxy);
        for (int y = y0; y <= y1; y++)
            for (int x = x0; x <= x1; x++) debut[y * g->nx + x + 1]++;
        total += (size_t)(x1 - x0 + 1) * (size_t)(y1 - y0 + 1);
    }
    if (total > g->cap_contenu) {
        g->contenu = xrealloc(g->contenu, total * sizeof *g->contenu);
        g->cap_contenu = total;
    }
    for (size_t c = 1; c <= nc; c++) debut[c] += debut[c - 1];

    /* debut[c] sert de curseur d'écriture : après le placement, il vaut
     * l'ancien debut[c+1], d'où le décalage final d'une case. */
    for (int i = 0; i < m->n; i++) {
        Rect b = boite_objet(&m->obj[i]);
        int x0 = cellule_x(g, b.minx), x1 = cellule_x(g, b.maxx);
        int y0 = cellule_y(g, b.miny), y1 = cellule_y(g, b.maxy);
        for (int y = y0; y <= y1; y++)
            for (int x = x0; x <= x1; x++) g->contenu[debut[y * g->nx + x]++] = i;
    }
    memmove(debut + 1, debut, nc * sizeof *debut);
    debut[0] = 0;
    g->nb_entrees = total;
}

static void *grille_construire(const Monde *m)
{
    Grille *g = xcalloc(1, sizeof *g);
    float c = taille_imposee;
    if (c <= 0) {
        double somme = 0;
        for (int i = 0; i < m->n; i++) somme += m->obj[i].r;
        float r_moyen = m->n > 0 ? (float)(somme / m->n) : 1.0f;
        c = sqrtf(m->largeur * m->hauteur / (float)(m->n > 0 ? m->n : 1));
        if (c < 2 * r_moyen) c = 2 * r_moyen;
    }
    g->taille = c;
    g->inv = 1.0f / c;
    g->nx = (int)ceilf(m->largeur / c);
    g->ny = (int)ceilf(m->hauteur / c);
    if (g->nx < 1) g->nx = 1;
    if (g->ny < 1) g->ny = 1;
    g->debut = xmalloc(((size_t)g->nx * g->ny + 1) * sizeof *g->debut);
    remplir(g, m);
    return g;
}

static void grille_liberer(void *s)
{
    Grille *g = s;
    free(g->debut);
    free(g->contenu);
    free(g);
}

/* Parcourt les cellules qui touchent la zone et garde les objets qui la
 * touchent (test exact du rectangle, ou du cercle si cercle != 0). */
static void parcourir(const Grille *g, const Monde *m, Rect zone,
                      int cercle, float cx, float cy, float rayon, ListeInt *res)
{
    int x0 = cellule_x(g, zone.minx), x1 = cellule_x(g, zone.maxx);
    int y0 = cellule_y(g, zone.miny), y1 = cellule_y(g, zone.maxy);
    for (int y = y0; y <= y1; y++) {
        for (int x = x0; x <= x1; x++) {
            int c = y * g->nx + x;
            for (int e = g->debut[c]; e < g->debut[c + 1]; e++) {
                int i = g->contenu[e];
                const Objet *o = &m->obj[i];
                Rect b = boite_objet(o);
                if (!rects_se_coupent(b, zone)) continue;
                /* Un gros objet est inscrit dans plusieurs cellules : on ne
                 * le compte que dans celle qui contient le coin bas-gauche
                 * de l'intersection objet ∩ zone (une seule cellule). */
                if (cellule_x(g, fmaxf(b.minx, zone.minx)) != x
                    || cellule_y(g, fmaxf(b.miny, zone.miny)) != y) continue;
                int dedans = cercle ? objet_dans_cercle(o, cx, cy, rayon)
                                    : objet_dans_rect(o, zone);
                if (dedans) liste_ajouter(res, i);
            }
        }
    }
}

static void grille_rect(void *s, const Monde *m, Rect zone, ListeInt *res)
{
    parcourir(s, m, zone, 0, 0, 0, 0, res);
}

static void grille_cercle(void *s, const Monde *m, float cx, float cy, float rayon, ListeInt *res)
{
    parcourir(s, m, boite_cercle(cx, cy, rayon), 1, cx, cy, rayon, res);
}

/* Propose au tas les objets dont le centre est dans la cellule (x, y). */
static void visiter(const Grille *g, const Monde *m, int x, int y,
                    float px, float py, TasK *tas)
{
    int c = y * g->nx + x;
    for (int e = g->debut[c]; e < g->debut[c + 1]; e++) {
        int i = g->contenu[e];
        const Objet *o = &m->obj[i];
        /* Chaque objet n'est examiné que dans la cellule de son centre. */
        if (cellule_x(g, o->x) != x || cellule_y(g, o->y) != y) continue;
        tas_proposer(tas, i, dist2(px, py, o->x, o->y));
    }
}

/* k plus proches voisins : on visite des anneaux de cellules de plus en
 * plus larges autour de la cellule du point, et on s'arrête dès que les
 * cellules restantes sont toutes plus loin que le k-ième candidat. */
static int grille_kppv(void *s, const Monde *m, float px, float py, int k, int *ids, float *d2)
{
    const Grille *g = s;
    if (k <= 0 || m->n == 0) return 0;
    TasK tas;
    tas_init(&tas, k, ids, d2);
    int qx = cellule_x(g, px), qy = cellule_y(g, py);

    for (int r = 0;; r++) {
        /* Anneau r : cellules à distance de Tchebychev exactement r. */
        for (int y = qy - r; y <= qy + r; y++) {
            if (y < 0 || y >= g->ny) continue;
            if (y == qy - r || y == qy + r) {          /* ligne complète */
                int x0 = qx - r < 0 ? 0 : qx - r;
                int x1 = qx + r >= g->nx ? g->nx - 1 : qx + r;
                for (int x = x0; x <= x1; x++) visiter(g, m, x, y, px, py, &tas);
            } else {                                    /* deux bords */
                if (qx - r >= 0) visiter(g, m, qx - r, y, px, py, &tas);
                if (qx + r < g->nx) visiter(g, m, qx + r, y, px, py, &tas);
            }
        }

        /* Distance minimale entre le point et les cellules pas encore
         * visitées : le côté le plus proche du carré déjà parcouru (les
         * côtés collés au bord du monde ne comptent pas). */
        float borne = INFINITY;
        int reste = 0;
        if (qx - r > 0)         { reste = 1; borne = fminf(borne, px - (float)(qx - r) * g->taille); }
        if (qx + r < g->nx - 1) { reste = 1; borne = fminf(borne, (float)(qx + r + 1) * g->taille - px); }
        if (qy - r > 0)         { reste = 1; borne = fminf(borne, py - (float)(qy - r) * g->taille); }
        if (qy + r < g->ny - 1) { reste = 1; borne = fminf(borne, (float)(qy + r + 1) * g->taille - py); }
        if (!reste) break;                              /* toute la grille est vue */
        /* Petite marge pour les arrondis du calcul des cellules. */
        float lim = borne - 0.001f * g->taille;
        if (tas.n == k && lim > 0 && tas_seuil(&tas) <= lim * lim) break;
    }
    return tas_finir(&tas);
}

/* Découpage pour le parallélisme : la part p traite les lignes de
 * cellules y = p, p + nb, p + 2 nb... */
static void grille_paires(void *s, const Monde *m, int part, int nb_parts, ListePaires *res)
{
    const Grille *g = s;
    for (int y = part; y < g->ny; y += nb_parts) {
        for (int x = 0; x < g->nx; x++) {
            int c = y * g->nx + x;
            int d = g->debut[c], f = g->debut[c + 1];
            for (int ea = d; ea < f; ea++) {
                int i = g->contenu[ea];
                Rect bi = boite_objet(&m->obj[i]);
                for (int eb = ea + 1; eb < f; eb++) {
                    int j = g->contenu[eb];
                    Rect bj = boite_objet(&m->obj[j]);
                    if (!rects_se_coupent(bi, bj)) continue;
                    /* Deux objets peuvent partager plusieurs cellules : la
                     * paire n'est comptée que dans la cellule du coin
                     * bas-gauche de l'intersection de leurs boîtes. */
                    if (cellule_x(g, fmaxf(bi.minx, bj.minx)) != x
                        || cellule_y(g, fmaxf(bi.miny, bj.miny)) != y) continue;
                    if (collision(&m->obj[i], &m->obj[j])) paires_ajouter(res, i, j);
                }
            }
        }
    }
}

static void grille_mettre_a_jour(void *s, const Monde *m)
{
    /* Le tri par dénombrement est en O(n + C) : reconstruire coûte moins
     * cher que de chercher quels objets ont changé de cellule. */
    remplir(s, m);
}

static size_t grille_memoire(const void *s)
{
    const Grille *g = s;
    return sizeof *g + ((size_t)g->nx * g->ny + 1) * sizeof *g->debut
         + g->cap_contenu * sizeof *g->contenu;
}

const Algo ALGO_GRILLE = {
    "grille", grille_construire, grille_liberer, grille_rect, grille_cercle,
    grille_kppv, grille_paires, grille_mettre_a_jour, grille_memoire
};
