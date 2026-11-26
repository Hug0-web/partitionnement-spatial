# Partitionnement spatial pour les jeux
 
Projet d'algorithmique en **C** : retrouver vite les objets proches dans un monde 2D qui contient des centaines à des dizaines de milliers d'objets, statiques ou en mouvement.
 
Le programme sait :
 
- trouver les objets présents dans une **zone** rectangulaire ou circulaire ;
- trouver les **k plus proches voisins** d'un point ;
- détecter **toutes les paires d'objets en collision** ;
- **mettre à jour** la structure quand les objets bougent.
Toutes les structures ont la même interface et sont comparées entre elles :
 
| Structure | Rôle dans le sujet |
|---|---|
| Force brute | Référence demandée |
| Grille uniforme | Grille demandée |
| Grille hachée (spatial hashing) | Variante de la grille, pour les mondes très étendus |
| kd-tree, en 3 variantes de mise à jour : réajusté, reconstruit, hybride | Structure arborescente demandée |
| Sweep and Prune | Bonus proposé par le sujet |
| Version parallèle (threads POSIX), pour toutes les structures | Bonus proposé par le sujet |
 
## Compilation et lancement
 
Il faut un compilateur C11 (clang ou gcc) et `make`. Python 3 et matplotlib ne servent qu'aux courbes, ils sont facultatifs.
 
```bash
cd src
make             # compile ./spatial
make test        # 404 vérifications : chaque structure contre la force brute
make svg         # figures des structures        -> ../figures/*.svg
make bench       # mesures principales            -> ../resultats/mesures.csv
make dynamique   # scène en mouvement             -> ../resultats/dynamique.csv, seuil_hybride.csv
make grille      # taille des cellules            -> ../resultats/grille.csv
make etendue     # monde de plus en plus étendu   -> ../resultats/etendue.csv
make parallele   # 1, 2, 4 et 8 fils d'exécution  -> ../resultats/parallele.csv
make tout        # tout ce qui précède (≈ 11 min)
make graphiques  # courbes PNG (matplotlib)       -> ../figures/*.png
make clean
```
 
Les mesures doivent tourner sur une machine éveillée et peu chargée. Sur un Mac, lancez-les avec `caffeinate -di make tout`, capot ouvert (voir section 4).
 
## Organisation
 
```
projet/
├── src/
│   ├── commun.c/.h     objets, géométrie, listes, tas des k voisins, hasard, chronomètre
│   ├── monde.c/.h      génération des scènes, déplacement des objets
│   ├── algos.c/.h      interface commune à toutes les structures
│   ├── naif.c          force brute
│   ├── grille.c/.h     grille uniforme
│   ├── hachage.c/.h    grille hachée (spatial hashing)
│   ├── kdtree.c/.h     kd-tree : réajusté, reconstruit ou hybride
│   ├── sap.c           Sweep and Prune
│   ├── parallele.c/.h  version parallèle (threads POSIX)
│   ├── bench.c/.h      mesures
│   ├── visu.c/.h       figures SVG
│   ├── tests.c         tests automatiques
│   ├── main.c
│   └── Makefile
├── graphiques.py       courbes à partir des CSV
├── resultats/          mesures (CSV, sortie console, durée de la campagne)
├── figures/            SVG des structures, PNG des mesures
├── rapport.pdf         ce document mis en page
└── README.md
```
 
Chaque structure fournit les mêmes fonctions : `construire`, `requete_rect`, `requete_cercle`, `kppv`, `paires`, `mettre_a_jour` et `memoire`. Les structures ne copient pas les objets : elles rangent des indices et relisent les positions au moment des requêtes. Elles ne sont que lues pendant les requêtes, ce qui permet de les interroger depuis plusieurs fils d'exécution à la fois, sans verrou.
 
## 1. Modèle
 
**Objets.** Chaque objet est un disque. Il a un centre, un rayon et une vitesse.
 
- Deux objets sont en collision si leurs disques se touchent.
- Un objet est dans une zone si son disque touche la zone.
- Pour les voisins, on mesure la distance jusqu'au centre de l'objet.
**Scènes.** Le monde est un carré de côté 10 √n. Plus il y a d'objets, plus le monde est grand. Il y a donc toujours à peu près un objet pour 100 unités². Comme ça, on voit l'effet du nombre d'objets sans que la scène devienne plus serrée.
 
On teste quatre scènes :
 
| | Tailles proches (rayon 1) | Tailles variées |
|---|---|---|
| **Uniforme** | objets placés au hasard partout | pareil |
| **Amas** | groupes de 250 objets qui couvrent 12 % du monde | pareil |
 
Avec des tailles variées, 95 % des objets ont un rayon entre 0,5 et 1,5. Les 5 % restants sont gros, avec un rayon entre 4 et 15.
 
On teste aussi :
 
- **une scène qui bouge** : à chaque pas, les objets avancent d'au plus 0,5 et rebondissent sur les bords ;
- **un très grand monde** : on agrandit le monde jusqu'à 32 fois en largeur. Les amas restent pareils, donc le monde est presque vide.
**Requêtes.** On choisit un objet au hasard et on regarde autour de lui, comme un personnage dans un jeu. On cherche les objets dans un carré de 20 × 20, dans un cercle de rayon 10, puis les 10 objets les plus proches. Chaque mesure se fait sur 1 000 requêtes.
 
## 2. Les algorithmes
 
Notations : n = nombre d'objets, m = nombre de résultats, C = nombre de cases, e = nombre de places occupées dans la grille (à peu près n si les objets sont petits), k = nombre de voisins.
 
| Opération | Force brute | Grille | Grille hachée | kd-tree | Sweep and Prune |
|---|---|---|---|---|---|
| Construction | O(1) | O(n + C + e) | O(n + e) | O(n log n) | O(n log n) |
| Zone | O(n) | O(cases touchées + m) | comme la grille | O(√n + m) | O(log n + tranche en x) |
| k voisins | O(n log k) | O(k) en moyenne | comme la grille | O(log n + k) en moyenne | O(√n) en moyenne* |
| Collisions | O(n²) | O(n) | comme la grille | O(n log n + paires) | O(n^1,5)* |
| Mise à jour | rien | on refait tout, O(n + C) | on refait tout, O(n + e) | O(n) ou O(n log n) | O(n + échanges) |
| Mémoire | 0 | 4 (C + 1) + 4e octets | 8 S + 4e octets (S = nombre de seaux, entre 2n et 4n) | 24 n octets | 8 n octets |
 
\* Une bande verticale du monde contient environ √n objets.
 
### Force brute ([naif.c](src/naif.c))
 
On compare chaque objet avec tous les autres. Pour les collisions, ça fait n(n − 1)/2 tests. C'est lent mais c'est sûr. On s'en sert pour vérifier les autres structures.
 
### Grille uniforme ([grille.c](src/grille.c))
 
On découpe le monde en cases carrées. On range chaque objet dans les cases qu'il touche. Un gros objet peut toucher plusieurs cases. Pour chercher autour d'un point, on regarde seulement les cases proches.
 
- **Stockage.** On utilise deux tableaux. `debut` dit où commence chaque case et `contenu` contient les objets. On remplit la grille en deux passes. D'abord on compte les objets de chaque case, puis on les range. Il n'y a pas de liste chaînée.
- **Pas de doublons.** Un objet rangé dans deux cases pourrait être trouvé deux fois. Pour l'éviter, on ne le garde que dans une case : celle qui contient le coin en bas à gauche de la partie commune entre l'objet et la zone. Pour une paire d'objets, on fait pareil avec la partie commune des deux objets. Ce coin est unique, donc il n'y a pas de doublon.
- **Voisins.** On regarde les cases autour du point, puis un anneau de cases plus large, et ainsi de suite. On s'arrête quand aucun objet plus loin ne peut être plus proche.
- **Quand les objets bougent.** On refait toute la grille. C'est rapide et plus simple que de déplacer les objets un par un.
- **Taille des cases.** Elle est choisie pour avoir environ un objet par case. Elle fait aussi au moins deux fois le rayon moyen. Voir la section 5.5.
### Grille hachée ([hachage.c](src/hachage.c))
 
C'est une grille sans limites. On n'a pas un tableau avec toutes les cases. À la place, chaque case est rangée dans un « seau » choisi par une formule (Teschner et al., 2003) :
 
h(x, y) = (73856093·x XOR 19349663·y) modulo S
 
S est le nombre de seaux. C'est la puissance de 2 juste au-dessus de 2n. La mémoire dépend donc du nombre d'objets et pas de la taille du monde.
 
- Deux cases peuvent tomber dans le même seau. Le test du coin en bas à gauche règle ce problème, car il écarte les objets des autres cases.
- Un objet peut toucher deux cases du même seau. Un tableau `dernier` évite de le ranger deux fois.
- Les cases font 10 fois le rayon moyen.
### kd-tree ([kdtree.c](src/kdtree.c))
 
C'est un arbre (Bentley, 1975). À chaque nœud, on coupe les objets en deux moitiés égales. On coupe au milieu (la médiane), dans la direction où les objets sont le plus étalés. On continue jusqu'à avoir au plus 8 objets par feuille.
 
- **Rectangles.** Chaque nœud garde un rectangle qui contient tous ses objets. Pendant une recherche, si le rectangle ne touche pas la zone, on ne descend pas dedans. Ça ressemble à un BVH.
- **Voisins.** On descend d'abord du côté le plus proche. On saute les parties de l'arbre qui sont trop loin.
- **Collisions.** Au lieu de faire une recherche par objet, on en fait une par feuille. Ça fait 8 fois moins de recherches. À 50 000 objets, ça prend 2,3 ms au lieu de 6,8 ms.
- **Quand les objets bougent**, j'ai essayé trois méthodes :
  - **réajustée** : on garde l'arbre et on recalcule juste les rectangles. C'est très rapide.
  - **reconstruite** : on refait tout l'arbre à chaque pas.
  - **hybride** : on réajuste, mais on reconstruit quand les rectangles deviennent trop grands. Ici, « trop grands » veut dire que la surface totale des feuilles a grossi de 50 %. Plus les rectangles sont grands, plus les recherches sont lentes. L'idée vient de MacDonald et Booth (1990). Le choix de 50 % est expliqué en section 5.4.
### Sweep and Prune ([sap.c](src/sap.c))
 
(Cohen et al., 1995.) On trie les objets de gauche à droite. Un objet ne peut toucher que les objets qui commencent avant qu'il se termine. On s'arrête donc dès qu'on tombe sur un objet qui commence plus loin.
 
- **Quand les objets bougent.** Ils bougent peu, donc la liste est presque déjà triée. Un tri par insertion la remet en ordre très vite.
- **Zones.** On cherche le début dans la liste triée, puis on parcourt une tranche.
- **Défaut.** Il ne regarde qu'un seul axe. Deux objets à la même position en x sont comparés même s'ils sont très loin en y.
### Version parallèle ([parallele.c](src/parallele.c))
 
On partage le travail entre plusieurs threads. Chaque thread prend une partie. Par exemple, avec P threads et la grille, le thread p prend les lignes p, p + P, p + 2P, etc. Chaque thread écrit dans sa propre liste et on les colle à la fin.
 
Les structures ne sont jamais modifiées pendant ce calcul. Il n'y a donc pas besoin de verrou.
 
## 3. Tests
 
`make test` lance 404 vérifications. Elles passent toutes.
 
- **Comparaison avec la force brute.** Sur les 4 scènes avec 3 000 objets, chaque structure doit trouver exactement les mêmes résultats que la force brute. On teste 400 rectangles, 400 cercles, 400 recherches de voisins et toutes les collisions. Ensuite on fait bouger les objets pendant 40 pas et on recommence.
- **Cas spéciaux.** Un monde vide, un seul objet, 200 objets au même endroit, un objet énorme qui couvre tout, des cases très petites ou très grandes.
- **Grille hachée difficile.** Seulement 16 seaux, ou un monde 8 fois plus large.
- **kd-tree hybride.** On vérifie qu'il se reconstruit bien tout seul (10 fois en 150 pas) et qu'il reste juste.
- **Version parallèle.** Avec 1, 3 et 8 threads, elle doit trouver la même chose que la version normale.
- **Les tests trouvent-ils vraiment les bugs ?** Pour le vérifier, j'ai mis exprès 8 bugs dans le code, un par un. Par exemple un doublon dans la grille ou une feuille oubliée dans le kd-tree. À chaque fois, entre 2 et 45 tests ont échoué.
- **Erreurs mémoire.** Aucune erreur avec AddressSanitizer et UndefinedBehaviorSanitizer.
Pour que la comparaison soit exacte, tous les algorithmes utilisent les mêmes fonctions de calcul ([commun.h](src/commun.h)). Je compile aussi avec `-ffp-contract=off`. Sans cette option, le compilateur peut faire les calculs un peu différemment d'un fichier à l'autre, et les arrondis ne sont plus les mêmes.
 
## 4. Comment j'ai mesuré
 
- **Machine.** Mac avec une puce M1 Pro et 16 Go de mémoire, compilé avec clang en `-O2`. Tout tourne sur un seul cœur, sauf la version parallèle.
- **Répétitions.** Chaque mesure dure au moins 10 ms. Je la fais 3 fois et je garde la plus rapide.
- **Données qui changent.** Si on refait le même calcul des milliers de fois, le processeur finit par le connaître par cœur. Il va alors jusqu'à 3,6 fois plus vite, ce qui n'est pas réaliste. J'utilise donc 8 versions de la scène l'une après l'autre, comme les images d'un jeu.
- **Mise en veille.** Ma première série de mesures était fausse, jusqu'à 100 fois trop lente. Le Mac se mettait en veille parce que le capot était fermé. Je l'ai refaite avec `caffeinate`, qui empêche la mise en veille. La durée est dans [resultats/temps.txt](resultats/temps.txt).
- **Précision.** Sans mouvement, les trois versions du kd-tree font exactement le même calcul. Leurs différences montrent donc l'erreur de mesure : environ 1,3 % en général et 13,5 % au pire.
- **Mémoire.** Je compte seulement la mémoire de la structure. Les objets prennent en plus 20 octets chacun.
## 5. Résultats
 
Tous les résultats sont dans [resultats/mesures.csv](resultats/mesures.csv) et [resultats/console.txt](resultats/console.txt). Les courbes sont dans `figures/*.png`.
 
### 5.1 Résultats à 50 000 objets
 
Pour le kd-tree, la mise à jour a trois valeurs : réajusté / reconstruit / hybride.
 
**Scène uniforme, tailles proches** (3 197 paires en collision, 5,8 objets par zone)
 
| | Construction | Zone rect. | Zone cercle | 10 voisins | Collisions | Mise à jour | Mémoire |
|---|---|---|---|---|---|---|---|
| Force brute | — | 110 µs | 96 µs | 82 µs | **2 784 ms** | — | 0 |
| Grille | 0,86 ms | **0,19 µs** | 0,24 µs | **0,91 µs** | **0,55 ms** | 0,81 ms | 477 Ko |
| Grille hachée | **0,84 ms** | 0,24 µs | **0,23 µs** | 1,13 µs | 0,69 ms | 0,83 ms | 1 305 Ko |
| kd-tree | 8,3 ms | 0,39 µs | 0,36 µs | 1,05 µs | 2,34 ms | **0,06** / 6,47 / 0,20 ms | 1 172 Ko |
| Sweep and Prune | 3,9 ms | 1,41 µs | 1,21 µs | 4,22 µs | 6,14 ms | 0,61 ms | **391 Ko** |
 
**Scène en amas, tailles variées** (125 793 paires, 52,6 objets par zone)
 
| | Construction | Zone rect. | Zone cercle | 10 voisins | Collisions | Mise à jour | Mémoire |
|---|---|---|---|---|---|---|---|
| Force brute | — | 111 µs | 95 µs | 82 µs | 2 783 ms | — | 0 |
| Grille | 0,95 ms | 1,72 µs | 1,75 µs | 1,81 µs | 5,97 ms | 0,91 ms | 553 Ko |
| Grille hachée | **0,80 ms** | 1,56 µs | 1,65 µs | 1,98 µs | **5,82 ms** | 0,76 ms | 1 321 Ko |
| kd-tree | 8,0 ms | **1,28 µs** | **1,15 µs** | **1,74 µs** | 7,99 ms | **0,06** / 5,70 / 0,19 ms | 1 172 Ko |
| Sweep and Prune | 4,0 ms | 4,62 µs | 3,83 µs | 4,44 µs | 13,23 ms | 0,63 ms | **391 Ko** |
 
La force brute met 2,8 s pour trouver les collisions. Dans un jeu à 60 images par seconde, on a 16,7 ms par image. La force brute est donc 167 fois trop lente. La grille fait la mise à jour et les collisions en 1,4 ms, soit 8 % du temps d'une image.
 
![Toutes les collisions selon n](figures/mesures_paires_ms.png)
 
### 5.2 Théorie et mesures
 
On passe de 5 000 à 50 000 objets, donc 10 fois plus. On compare ce que prévoit la complexité avec ce qu'on mesure (scène uniforme, tailles proches).
 
| Opération | Complexité | Prévu | Mesuré |
|---|---|---|---|
| Force brute, collisions | O(n²) | × 100 | × 101 (27,5 ms → 2 784 ms) |
| Force brute, zone | O(n) | × 10 | × 10,0 |
| Grille, construction | O(n) | × 10 | × 10,9 |
| Grille, collisions | O(n) | × 10 | × 12,3 |
| Grille, zone | O(1) | × 1 | × 1,10 |
| Grille hachée, construction | O(n) | × 10 | × 10,8 |
| kd-tree, construction | O(n log n) | × 12,7 | × 12,1 |
| kd-tree, zone | O(log n) pour une petite zone | × 1,27 | × 1,45 |
| kd-tree, collisions | O(n log n) | × 12,7 | × 10,2 |
| kd-tree, réajustement | O(n) | × 10 | × 10,9 |
| Sweep and Prune, construction | O(n log n) | × 12,7 | × 12,8 |
| Sweep and Prune, zone | O(√n) | × 3,2 | × 2,8 |
| Sweep and Prune, collisions | O(n^1,5) | × 31,6 | × 27,1 |
 
Les mesures sont proches de la théorie. Quand c'est un peu plus lent que prévu, c'est à cause du cache. À 50 000 objets, les données ne tiennent plus dans la mémoire rapide du processeur.
 
Pour le Sweep and Prune, le temps monte plus vite que n. C'est normal : plus il y a d'objets, plus chaque objet a de voisins sur l'axe x.
 
### 5.3 Effet des amas et des gros objets
 
Temps pour trouver toutes les collisions à 50 000 objets (ms) :
 
| Scène | Paires | Grille | Grille hachée | kd-tree | Sweep and Prune |
|---|---|---|---|---|---|
| Uniforme, tailles proches | 3 197 | **0,55** | 0,69 | 2,34 | 6,14 |
| Uniforme, tailles variées | 12 574 | **1,03** | 1,07 | 3,27 | 8,81 |
| Amas, tailles proches | 31 569 | **1,83** | 1,95 | 4,09 | 8,98 |
| Amas, tailles variées | 125 793 | 5,97 | **5,82** | 7,99 | 13,23 |
 
- **Les grilles sont les plus rapides partout pour les collisions.** Mais elles n'aiment pas les amas : elles passent de 0,55 à 6,0 ms, environ 11 fois plus lent. Dans un amas, chaque case contient beaucoup d'objets, et il faut tous les comparer entre eux.
- **Le kd-tree supporte mieux les amas.** Il coupe toujours au milieu, donc chaque feuille garde 8 objets. Il ne ralentit que 3,4 fois alors qu'il y a 39 fois plus de paires.
- **Dans les amas, le kd-tree est le meilleur pour les zones et les voisins.** Par exemple, pour une zone : 1,28 µs contre 1,72 µs pour la grille. Sans amas, c'est la grille qui gagne.
- **Les gros objets gênent le Sweep and Prune.** La tranche à parcourir devient plus large. Sa zone passe de 1,41 à 2,92 µs.
- **Ils gênent aussi un peu la grille**, car un gros objet est rangé dans plusieurs cases.
### 5.4 Objets en mouvement
 
10 000 objets et 300 pas. À chaque pas, on déplace les objets, on met à jour la structure et on cherche les collisions. La force brute n'est pas dans le tableau car elle prend 110 ms par pas. Les données sont dans [resultats/dynamique.csv](resultats/dynamique.csv).
 
| Uniforme, tailles proches | Mise à jour | Collisions | Collisions pas 1-10 | Collisions pas 291-300 | **Total par pas** | Reconstructions |
|---|---|---|---|---|---|---|
| Grille | 0,16 ms | 0,11 ms | 0,11 ms | 0,11 ms | **0,27 ms** | — |
| Grille hachée | 0,16 ms | 0,13 ms | 0,12 ms | 0,13 ms | 0,29 ms | — |
| kd-tree réajusté | 0,01 ms | 5,23 ms | 0,50 ms | **12,85 ms** | 5,24 ms | 0 |
| kd-tree reconstruit | 1,21 ms | 0,48 ms | 0,50 ms | 0,49 ms | 1,69 ms | 300 |
| kd-tree hybride | 0,06 ms | 0,57 ms | 0,51 ms | 0,62 ms | **0,63 ms** | 11 |
| Sweep and Prune | 0,12 ms | 0,61 ms | 0,62 ms | 0,60 ms | 0,73 ms | — |
 
| Amas, tailles variées | Mise à jour | Collisions | Collisions pas 1-10 | Collisions pas 291-300 | **Total par pas** | Reconstructions |
|---|---|---|---|---|---|---|
| Grille | 0,18 ms | 0,50 ms | 1,02 ms | 0,29 ms | 0,68 ms | — |
| Grille hachée | 0,15 ms | 0,49 ms | 1,01 ms | 0,27 ms | **0,64 ms** | — |
| kd-tree réajusté | 0,01 ms | 10,16 ms | 1,63 ms | **21,53 ms** | 10,18 ms | 0 |
| kd-tree reconstruit | 1,15 ms | 1,02 ms | 1,54 ms | 0,76 ms | 2,17 ms | 300 |
| kd-tree hybride | 0,06 ms | 1,14 ms | 1,61 ms | 0,96 ms | 1,20 ms | 11 |
| Sweep and Prune | 0,12 ms | 1,22 ms | 1,54 ms | 1,05 ms | 1,34 ms | — |
 
- **Réajuster seulement, c'est un piège.** C'est très rapide au début. Mais les objets d'une même feuille s'éloignent les uns des autres et les rectangles grossissent. Après 300 pas, les collisions sont 25 fois plus lentes (0,50 → 12,85 ms). On le voit bien sur les figures plus bas.
- **L'hybride règle ce problème.** Il ne reconstruit que 11 fois en 300 pas. Il est 2,7 fois plus rapide que reconstruire à chaque pas et 8 fois plus rapide que réajuster seulement.
- **Le Sweep and Prune profite de sa liste presque triée.** Le tri lui prend 0,12 ms au lieu de 0,68 ms pour un tri complet.
- **Les grilles restent les plus rapides**, car les refaire ne coûte que 0,16 ms.
- **En amas, les collisions baissent avec le temps**, car les groupes se dispersent.
![Coût par pas en scène dynamique](figures/dynamique.png)
 
**Seuil de l'hybride.** J'ai essayé plusieurs seuils ([resultats/seuil_hybride.csv](resultats/seuil_hybride.csv)). Les valeurs sont en ms par pas, avec le nombre de reconstructions entre parenthèses.
 
| Seuil | 1,05 | 1,1 | 1,25 | **1,5** | 2 | 3 | 5 |
|---|---|---|---|---|---|---|---|
| Uniforme, tailles proches | 0,66 (37) | 0,64 (27) | **0,62** (16) | 0,65 (11) | 0,70 (8) | 0,79 (5) | 1,01 (3) |
| Amas, tailles variées | 1,17 (36) | **1,13** (25) | 1,14 (15) | 1,20 (11) | 1,31 (7) | 1,54 (5) | 1,94 (3) |
 
Si le seuil est trop bas, on reconstruit trop souvent. S'il est trop haut, l'arbre devient mauvais. Entre 1,1 et 1,5, ça change peu. J'ai gardé 1,5 : c'est à moins de 7 % du meilleur et ça ne fait que 11 reconstructions.
 
![kd-tree réajusté après 300 pas](figures/kdtree_reajuste_pas300.svg)
 
![kd-tree reconstruit au même instant](figures/kdtree_reconstruit_pas300.svg)
 
### 5.5 Taille des cases de la grille
 
10 000 objets, temps des collisions en ms. Le choix automatique donne des cases de 10 ([resultats/grille.csv](resultats/grille.csv)).
 
| Taille des cases | 1 | 2 | 5 | 10 (auto) | 15 | 20 | 30 | 50 |
|---|---|---|---|---|---|---|---|---|
| Uniforme, tailles proches | 0,89 | 0,36 | 0,19 | 0,10 | 0,076 | **0,075** | 0,10 | 0,21 |
| Uniforme, tailles variées | 1,53 | 0,67 | 0,32 | 0,20 | 0,16 | **0,15** | 0,17 | 0,28 |
| Amas, tailles variées | 3,80 | 2,27 | 1,70 | 1,41 | **1,34** | 1,41 | 1,64 | 2,08 |
| Mémoire (uniforme, tailles proches) | 4,2 Mo | 1,1 Mo | 233 Ko | 96 Ko | 68 Ko | 57 Ko | 49 Ko | 44 Ko |
 
- **Cases trop petites.** Un objet est rangé dans beaucoup de cases. C'est 9 fois plus lent et ça prend 44 fois plus de mémoire.
- **Cases trop grandes.** Il y a trop d'objets dans chaque case. On revient presque à la force brute.
- **Le meilleur choix est entre 15 et 20.** Le choix automatique (10) est 35 % plus lent en scène uniforme et environ 5 % plus lent en amas. Il reste correct, mais pour un vrai jeu il vaut mieux tester.
### 5.6 Très grand monde : grille ou grille hachée ?
 
10 000 objets en amas avec des tailles variées. On agrandit le monde jusqu'à 32 fois en largeur, donc 1 024 fois en surface ([resultats/etendue.csv](resultats/etendue.csv)).
 
| Agrandissement du monde | 1 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| Grille, taille auto : mémoire | 110 Ko | 85 Ko | 81 Ko | 80 Ko | 79 Ko |
| Grille, taille auto : collisions | 1,12 ms | 1,23 ms | 1,59 ms | 1,81 ms | 1,91 ms |
| Grille, cases de 15 : mémoire | 76 Ko | 337 Ko | 1,2 Mo | 4,5 Mo | **17,8 Mo** |
| Grille, cases de 15 : collisions | 1,06 ms | 0,94 ms | 1,06 ms | 1,59 ms | 3,75 ms |
| **Grille hachée : mémoire** | 315 Ko | 316 Ko | 315 Ko | 316 Ko | **315 Ko** |
| **Grille hachée : collisions** | 1,04 ms | 0,90 ms | 0,88 ms | 0,90 ms | **0,87 ms** |
| kd-tree : mémoire / collisions | 235 Ko / 1,61 ms | 235 Ko / 1,40 ms | 235 Ko / 1,34 ms | 235 Ko / 1,34 ms | 235 Ko / 1,35 ms |
 
Avec une grille normale, on a deux choix et aucun n'est bon :
 
- **Petites cases.** Il faut beaucoup de cases pour couvrir le monde. La mémoire passe de 76 Ko à 17,8 Mo. Les collisions sont 3,5 fois plus lentes, car il faut parcourir plein de cases vides.
- **Grandes cases.** Chaque amas tient dans une ou deux cases. Les collisions sont 1,7 fois plus lentes.
**La grille hachée n'a pas ce problème.** Sa mémoire et son temps ne changent pas. C'est la plus rapide dans un grand monde. Le kd-tree ne change pas non plus, mais il est 1,5 fois plus lent.
 
Par contre, dans un monde normal, la grille hachée prend 2,7 fois plus de mémoire que la grille (1,3 Mo contre 0,48 Mo à 50 000 objets).
 
![Monde étendu](figures/etendue.png)
 
### 5.7 Version parallèle
 
50 000 objets, de 1 à 8 threads ([resultats/parallele.csv](resultats/parallele.csv)). Temps pour trouver toutes les collisions :
 
| | Uniforme : 1 → 8 threads | Gain | Amas : 1 → 8 threads | Gain |
|---|---|---|---|---|
| Force brute | 2 776 → 373 ms | **× 7,4** | 2 794 → 385 ms | **× 7,3** |
| Grille | 0,55 → 0,13 ms | × 4,2 | 5,68 → 0,88 ms | × 6,5 |
| Grille hachée | 0,68 → 0,18 ms | × 3,8 | 5,48 → 0,91 ms | × 6,0 |
| kd-tree | 2,29 → 0,44 ms | × 5,2 | 7,73 → 1,26 ms | × 6,1 |
| Sweep and Prune | 6,33 → 1,16 ms | × 5,5 | 13,88 → 2,14 ms | × 6,5 |
 
Quand chaque objet cherche ses 10 voisins (50 000 recherches) : la grille passe de 47 à 11 ms, la grille hachée de 58 à 13 ms et le kd-tree de 56 à 17 ms.
 
- **Plus il y a de travail, plus le parallélisme aide.** La force brute va presque 8 fois plus vite. Pour la grille en scène uniforme, le calcul est déjà très court (0,55 ms). Créer les threads prend alors une bonne partie du temps, et on ne gagne que × 4,2.
- **Exemple pour une IA.** Si chaque objet cherche ses 10 voisins, il faut 47 ms sur un cœur. C'est trop pour une seule image. Avec 8 cœurs on descend à 11 ms, et ça rentre.
![Accélération des collisions](figures/parallele.png)
 
### 5.8 Figures
 
`make svg` crée 8 images qu'on peut ouvrir dans un navigateur.
 
| Fichier | Ce qu'on voit |
|---|---|
| `figures/grille_uniforme.svg`, `figures/grille_amas.svg` | La grille. Les gros objets touchent plusieurs cases. |
| `figures/kdtree_uniforme.svg`, `figures/kdtree_amas.svg` | Les coupes du kd-tree et les rectangles des feuilles en orange. Les coupes sont plus serrées dans les amas. |
| `figures/kdtree_reajuste_pas000.svg`, `…pas100.svg`, `…pas300.svg` | Le kd-tree réajusté après 0, 100 et 300 pas. Les rectangles grossissent beaucoup. |
| `figures/kdtree_reconstruit_pas300.svg` | Le kd-tree reconstruit au même moment. Les rectangles restent petits. |
 
Couleurs : rouge pour les objets en collision, bleu pour les objets dans la zone (en pointillés), vert pour les 8 voisins les plus proches du point vert. La grille et le kd-tree trouvent les mêmes résultats.
 
![kd-tree sur la scène en amas](figures/kdtree_amas.svg)
 
## 6. Conclusion : quel algorithme pour quel jeu ?
 
| Situation | Meilleur choix | Pourquoi |
|---|---|---|
| Beaucoup de petits objets qui bougent : shoot'em up, particules, foules | **Grille** | La plus rapide. On la refait à chaque image en 0,16 ms. Il faut bien choisir la taille des cases. |
| Monde ouvert très grand et presque vide | **Grille hachée** | Sa mémoire ne grossit pas avec le monde. C'est la plus rapide dans ce cas. |
| Décor qui ne bouge pas, objets en groupes ou de tailles très différentes | **kd-tree** | Il s'adapte aux zones pleines et vides. On ne le construit qu'une fois. |
| Un arbre avec des objets qui bougent | **kd-tree hybride** | 2,7 fois plus rapide que tout reconstruire, et bien mieux que réajuster seulement. |
| Objets qui bougent surtout dans une direction : jeu de plateforme, course | **Sweep and Prune** | Mise à jour très rapide et peu de mémoire. À éviter si beaucoup d'objets ont la même position en x. |
| Beaucoup d'IA qui cherchent leurs voisins | **Grille + parallèle** | 47 ms sur un cœur, 11 ms sur 8 cœurs. |
| Quelques centaines d'objets | Force brute ou grille | 0,28 ms pour 500 objets. Mais déjà 110 ms pour 10 000. |
 
En résumé :
 
- **Quand tout bouge**, la grille gagne parce qu'elle est simple.
- **Quand rien ne bouge ou que les objets sont en groupes**, le kd-tree gagne.
- **À éviter** : la force brute avec des milliers d'objets.
## 7. Limites
 
- Tout est en 2D avec des disques. En 3D le principe reste le même.
- Pour les voisins, on mesure jusqu'au centre de l'objet et pas jusqu'à son bord.
- Le kd-tree hybride reste simple. Les vrais moteurs physiques utilisent des arbres plus poussés (BVH dynamiques).
- Le Sweep and Prune n'utilise qu'un seul axe. On pourrait en utiliser plusieurs.
- En parallèle, le travail est partagé de façon fixe. Dans les amas, certains threads ont plus de travail que les autres. La construction des structures n'est pas parallèle (voir Karras, 2012).
- Tout a été mesuré sur une seule machine. Les temps seront différents ailleurs, mais pas les conclusions.
## Références
 
- J. L. Bentley (1975), *Multidimensional binary search trees used for associative searching*, Communications of the ACM. L'origine du kd-tree.
- R. A. Finkel et J. L. Bentley (1974), *Quad trees: a data structure for retrieval on composite keys*, Acta Informatica.
- J. D. Cohen, M. C. Lin, D. Manocha, M. Ponamgi (1995), *I-COLLIDE: an interactive and exact collision detection system for large-scale environments*. Sweep and Prune.
- M. Teschner et al. (2003), *Optimized Spatial Hashing for Collision Detection of Deformable Objects*. Grille hachée et fonction de hachage utilisée dans [hachage.c](src/hachage.c).
- J. D. MacDonald et K. S. Booth (1990), *Heuristics for ray tracing using space subdivision*, The Visual Computer. Idée utilisée pour le kd-tree hybride.
- T. Karras (2012), *Maximizing Parallelism in the Construction of BVHs, Octrees, and k-d Trees*.
- G. Marsaglia (2003), *Xorshift RNGs*, Journal of Statistical Software. Générateur aléatoire utilisé dans [commun.c](src/commun.c).
 
