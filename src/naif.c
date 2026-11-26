/*
 * Force brute : aucune structure, on teste tous les objets.
 * Sert de référence pour les tests et les mesures.
 *
 *   Construction   O(1)
 *   Zone           O(n)
 *   k voisins      O(n log k)
 *   Collisions     O(n²)  : n(n-1)/2 tests
 *   Mise à jour    O(1)   : rien à faire
 *   Mémoire        O(1)
 */
#include "algos.h"

#include <stdlib.h>

typedef struct { int vide; } Naif;

static void *naif_construire(const Monde *m)
{
    (void)m;
    return xcalloc(1, sizeof(Naif));
}

static void naif_liberer(void *s)
{
    free(s);
}

static void naif_rect(void *s, const Monde *m, Rect zone, ListeInt *res)
{
    (void)s;
    for (int i = 0; i < m->n; i++) {
        if (objet_dans_rect(&m->obj[i], zone)) liste_ajouter(res, i);
    }
}

static void naif_cercle(void *s, const Monde *m, float cx, float cy, float rayon, ListeInt *res)
{
    (void)s;
    for (int i = 0; i < m->n; i++) {
        if (objet_dans_cercle(&m->obj[i], cx, cy, rayon)) liste_ajouter(res, i);
    }
}

static int naif_kppv(void *s, const Monde *m, float x, float y, int k, int *ids, float *d2)
{
    (void)s;
    if (k <= 0) return 0;
    TasK tas;
    tas_init(&tas, k, ids, d2);
    for (int i = 0; i < m->n; i++) {
        tas_proposer(&tas, i, dist2(x, y, m->obj[i].x, m->obj[i].y));
    }
    return tas_finir(&tas);
}

/* Découpage pour le parallélisme : la part p traite i = p, p + nb, p + 2 nb...
 * (en alternance, car les premiers i ont plus de j à tester). */
static void naif_paires(void *s, const Monde *m, int part, int nb_parts, ListePaires *res)
{
    (void)s;
    for (int i = part; i < m->n; i += nb_parts) {
        for (int j = i + 1; j < m->n; j++) {
            if (collision(&m->obj[i], &m->obj[j])) paires_ajouter(res, i, j);
        }
    }
}

static void naif_mettre_a_jour(void *s, const Monde *m)
{
    (void)s;
    (void)m;
}

static size_t naif_memoire(const void *s)
{
    (void)s;
    return 0;
}

const Algo ALGO_NAIF = {
    "naif", naif_construire, naif_liberer, naif_rect, naif_cercle,
    naif_kppv, naif_paires, naif_mettre_a_jour, naif_memoire
};
