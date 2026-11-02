# Courbes par splines cubiques : interpolation, racines et géométrie

## Présentation du projet

Ce projet construit une courbe lisse qui passe par des points donnés (une **spline cubique naturelle**), puis en étudie les propriétés avec des outils de calcul numérique : résolution de systèmes linéaires, polynômes, recherche de racines et géométrie du plan. On peut par exemple demander où la courbe coupe une droite, quel est son point le plus proche d'un point donné, quelle est sa longueur ou l'aire sous la courbe. Inspiré du tracé de routes et de trajectoires dans les jeux vidéo, il montre aussi pourquoi un seul polynôme de grand degré est une mauvaise idée (phénomène de Runge).

Tout le cœur du calcul est fait à la main, sans bibliothèque numérique (pas de numpy). `matplotlib` ne sert qu'à l'affichage et aux graphiques.



## Organisation du travail

Projet réalisé seul. Pour garder la même structure que les projets à trois, le travail est découpé en trois parties indépendantes.

**Partie commune (socle) :** choix du sujet, algèbre linéaire à la main (`algebre.py`), polynômes (`polynomes.py`), interpolation et spline (`interpolation.py`), affichage interactif à la souris (`main.py`).

**Approfondissement en trois parties :**

- **Partie 1 : conditionnement et pivot** (`mesure_conditionnement.py`, « Résultat 1 »). Étude de l'erreur de Gauss sur un système de Vandermonde, avec et sans pivot, contre une référence exacte en `Fraction`.
- **Partie 2 : Runge et coût de calcul** (`mesure_runge.py`, « Résultat 2 »). Comparaison polynôme à noeuds équidistants, polynôme à noeuds de Tchebychev et spline ; comparaison des temps entre Gauss et l'algorithme de Thomas.
- **Partie 3 : calcul exact et cas limites** (`verification_exacte.py`, « Résultat 3 »). Bissection contre Newton, polynôme de Wilkinson, racines multiples, comptage exact des racines par suites de Sturm, vérification des intersections courbe/droite.

**Tests :** `test_courbes.py` (26 tests) et les trois scripts de mesure ont été exécutés.

**TP :** les solutions des TP sont dans `tp/`.

## Structure du dépôt

```
tp/                          solutions des TP
projet/
  main.py                    affichage interactif (souris)
  algebre.py                 Gauss (avec/sans pivot), Thomas, inverse, conditionnement
  polynomes.py               Horner, dérivée, division, bissection, Newton, racines, Sturm
  interpolation.py           Vandermonde, différences divisées, Tchebychev, spline cubique
  geometrie.py               intersection courbe/droite, point le plus proche, longueur, aire
  mesure_conditionnement.py  partie 1
  mesure_runge.py            partie 2
  verification_exacte.py     partie 3
  test_courbes.py            tests automatiques
  resultats/                 tableaux (.csv) et graphiques (.png) produits par les scripts
```

## Comment lancer

Python 3 et matplotlib suffisent : `pip install matplotlib`. Depuis le dossier `projet/` :

| Commande | Effet |
|---|---|
| `python main.py` | ouvre la fenêtre interactive |
| `python mesure_conditionnement.py` | étude du pivot et du conditionnement (tableau + graphique) |
| `python mesure_runge.py` | Runge et temps de calcul (tableaux + graphiques) |
| `python verification_exacte.py` | racines, calcul exact et cas limites |
| `python -m unittest test_courbes -v` | lance les 26 tests |

**Ce qu'il faut observer avec `main.py` :** la courbe bleue passe par les carrés noirs. Une droite horizontale grise suit la souris ; ses intersections avec la courbe sont des points rouges. Le triangle vert est le point de la courbe le plus proche de la souris, relié à celle-ci par un trait vert. Le titre affiche la longueur de la courbe, l'aire sous la courbe, le nombre d'intersections et la distance à la courbe. Clic gauche : ajoute un point (refusé si son abscisse est trop proche d'une autre). Clic droit : supprime le point le plus proche (il en reste au moins 3). Touche `r` : revient aux points de départ.

## Notions du module mobilisées

- **Flottants et calcul exact :** erreurs d'arrondi dans Gauss, évaluation de Horner, tolérances, comparaison avec `Fraction` (calcul exact).
- **Arithmétique :** calculs exacts sur les rationnels, division euclidienne de polynômes, suites de Sturm.
- **Polynômes et racines :** schéma de Horner, dérivée, bissection, Newton, isolation des racines par les points critiques.
- **Matrices et systèmes linéaires :** élimination de Gauss avec pivot partiel, matrice de Vandermonde, conditionnement, système tridiagonal de la spline (algorithme de Thomas).
- **Géométrie du plan :** intersection courbe/droite (racines d'un polynôme de degré 3), plus court chemin d'un point à la courbe (racines d'un polynôme de degré 5), longueur d'arc, aire sous la courbe.
- **Mesure du coût :** temps d'exécution selon le nombre de points.

## Méthode

Sur chaque intervalle `[x_i, x_{i+1}]`, la spline est un polynôme de degré 3 en `t = x - x_i`. On impose que la courbe passe par les points, que la pente et la courbure soient continues, et que la courbure soit nulle aux deux bouts (spline « naturelle »). Les courbures inconnues vérifient un **système tridiagonal**, résolu soit par Thomas (en O(n)), soit par Gauss (même système en matrice pleine).

Pour **couper la courbe avec une droite** `y = m x + p`, on forme sur chaque segment le polynôme `S_i(t) - (m (x_i + t) + p)` et on cherche ses racines dans `[0, largeur]`. Pour trouver le **point le plus proche** d'un point `(p_x, p_y)`, la dérivée de la distance au carré est un polynôme de degré 5 en `t` ; on teste ses racines et les deux extrémités de chaque segment.

Les racines d'un polynôme sont isolées ainsi : entre deux points critiques consécutifs (racines de la dérivée, trouvées récursivement), le polynôme est monotone, donc il a au plus une racine, trouvée par bissection. Un point critique où le polynôme est quasi nul est pris comme racine double.

## Résultats

Les nombres ci-dessous viennent des scripts du dépôt (aucun aléa hors graine fixée). Les temps dépendent de la machine ; ceux du tableau 2b sont donnés à titre d'exemple.

### Résultat 1 : conditionnement et pivot (`mesure_conditionnement.py`)

On interpole `f(x) = 1/(1+x²)` sur `n` points équidistants de `[0, 1]` en résolvant le système de Vandermonde. L'erreur est l'erreur relative maximale sur les coefficients par rapport à la solution exacte calculée avec `Fraction`.

| n | cond(V) | erreur (pivot) | erreur (sans pivot) |
|---|---|---|---|
| 5 | 1,7e3 | 4,4e-15 | 4,4e-15 |
| 8 | 8,1e5 | 3,7e-13 | 1,1e-12 |
| 10 | 4,8e7 | 1,2e-10 | 1,6e-11 |
| 12 | 3,2e9 | 1,1e-9 | 1,7e-9 |
| 14 | 2,0e11 | 7,8e-8 | 7,0e-8 |
| 16 | 1,2e13 | 6,7e-6 | 7,1e-7 |
| 18 | 8,1e14 | 2,7e-4 | 1,7e-5 |
| 20 | 4,9e16 | 2,7e-2 | 6,0e-3 |

Cas d'un petit pivot : `[[eps, 1], [1, 1]] x = [1, 2]`, erreur sur `x0` par rapport à la solution exacte.

| eps | erreur avec pivot | erreur sans pivot |
|---|---|---|
| 1e-5 | 0 | 3,8e-12 |
| 1e-10 | 0 | 8,3e-8 |
| 1e-15 | 2,2e-16 | 8,0e-4 |
| 1e-20 | 0 | 1,0 |

**Interprétation.** Le conditionnement de la matrice de Vandermonde explose (il est multiplié par 8 environ à chaque point ajouté) et l'erreur le suit : à `n = 20`, `cond(V)` est de l'ordre de `1/ε_machine` et on n'a presque plus aucun chiffre correct. Ce n'est pas un défaut de l'algorithme : c'est le problème lui-même qui est mal posé, et le pivot ne peut pas le réparer (sur Vandermonde, avec ou sans pivot, l'erreur est du même ordre). Le pivot sert à autre chose : éviter de diviser par un très petit pivot, comme le montre le système 2×2, où l'absence de pivot donne une réponse entièrement fausse pour `eps = 1e-20`. Conclusion : on évite de passer par les coefficients de Vandermonde et on utilise la forme de Newton ou une spline.

### Résultat 2 : Runge et coût de calcul (`mesure_runge.py`)

**2a. Précision.** Erreur maximale (sur 2001 points) pour `f(x) = 1/(1+25x²)` sur `[-1, 1]`.

| n | polynôme équidistant | polynôme Tchebychev | spline |
|---|---|---|---|
| 5 | 4,4e-1 | 4,0e-1 | 2,8e-1 |
| 9 | 1,0 | 1,7e-1 | 5,6e-2 |
| 13 | 3,7 | 6,9e-2 | 6,9e-3 |
| 17 | 14,4 | 3,3e-2 | 3,7e-3 |
| 21 | 59,8 | 1,5e-2 | 3,2e-3 |
| 31 | 2,4e3 | 2,1e-3 | 8,2e-4 |
| 41 | 1,0e5 | 2,9e-4 | 2,8e-4 |

**2b. Temps de construction de la spline** (en ms, exemple de résultat) :

| Points | Thomas | Gauss |
|---|---|---|
| 10 | 0,009 | 0,025 |
| 20 | 0,017 | 0,072 |
| 40 | 0,027 | 0,232 |
| 80 | 0,061 | 0,825 |
| 160 | 0,112 | 3,143 |

**Interprétation.** Avec des noeuds équidistants, l'erreur du polynôme **augmente** avec `n` (elle dépasse 1 dès `n = 9`, soit plus que la valeur maximale de la fonction) : c'est le phénomène de Runge. Les noeuds de Tchebychev et la spline corrigent le problème, la spline étant la plus précise pour les petits `n`. Son erreur n'est pas strictement décroissante (`n = 15` puis `n = 17`) : elle dépend de la position des noeuds par rapport au pic de la fonction. Côté coût, quand on double le nombre de points, Thomas double son temps (coût linéaire) alors que Gauss le multiplie par environ 3,6 (coût proche du carré, car notre élimination saute les facteurs nuls, sinon on aurait un coût cubique). À 160 points, Thomas est environ 28 fois plus rapide.

### Résultat 3 : racines, calcul exact et cas limites (`verification_exacte.py`)

**3a. Bissection contre Newton** (tolérance 1e-12) :

| Cas | Bissection | Newton |
|---|---|---|
| x² - 2 sur [1, 2] | 40 it. | 5 it. |
| x³ - 2x - 5 sur [2, 3] | 40 it. | 6 it. |
| (x-1)²(x+2), racine double | impossible (pas de changement de signe) | 30 it. |
| x³ - 2x + 2 depuis x0 = 0 | 42 it. | échec : cycle 0 → 1 → 0 |

**3b. Polynôme de Wilkinson** `(x-1)(x-2)...(x-n)`, coefficients arrondis en flottants :

| n | racines exactes | trouvées | erreur max |
|---|---|---|---|
| 5 | 5 | 5 | 8,8e-13 |
| 10 | 10 | 10 | 1,4e-10 |
| 15 | 15 | 15 | 8,3e-7 |
| 17 | 17 | 17 | 2,7e-5 |
| 20 | 20 | 20 | 9,9e-3 |

**3c. Cas limites** dans `[0, 5]` (nombre de racines distinctes, exact par Sturm contre flottants) :

| Polynôme | Exact | Flottants |
|---|---|---|
| (x-1)(x-3) | 2 | 2 |
| (x-1)²(x-3), racine double | 2 | 2 |
| (x-1)³, racine triple | 1 | 1 |
| (x-1)(x-1,000001) | 2 | 2 |
| (x-1)(x-1-1e-12) | 2 | **1** |
| x² + 1 | 0 | 0 |

**3d. Intersections courbe/droite.** Spline de 8 points (ordonnées entières tirées au hasard, graine 2024) : sur 46 niveaux horizontaux, dont 6 situés exactement sur un noeud, le nombre d'intersections trouvé en flottants est égal au nombre exact (Sturm sur `Fraction`) dans tous les cas.

**Interprétation.** Newton est bien plus rapide quand il converge, mais il peut échouer (cycle) ou ralentir (racine double : convergence seulement linéaire). La bissection est lente mais sûre dès qu'il y a un changement de signe. Pour le polynôme de Wilkinson, une toute petite erreur d'arrondi sur les coefficients (de l'ordre de 1e-16 en relatif) déplace les racines jusqu'à 1e-2 à `n = 20` : le problème est très sensible, ce n'est pas l'algorithme qui est en cause. Enfin, deux racines distantes de 1e-12 sont confondues en flottants, mais les suites de Sturm en calcul exact les distinguent.

## Limites connues

- La courbe est le graphe d'une fonction : les abscisses des points doivent être strictement croissantes (pas de boucles ni de courbe fermée).
- Spline naturelle uniquement (courbure nulle aux extrémités) ; pas de spline périodique ni de contraintes sur les pentes.
- Les tolérances sont fixes : 1e-12 pour les pivots et la bissection, 1e-9 pour fusionner des racines proches ou détecter une racine double. Deux racines plus proches que cela sont confondues (3c).
- Le calcul exact (Sturm) ne marche qu'avec des `Fraction` ; les coefficients de la spline en flottants ne sont eux-mêmes déjà plus exacts.
- Le coût de `racines_reelles` est élevé pour des polynômes de grand degré (récursion sur les dérivées) ; il reste raisonnable ici (degrés 3 à 5 pour la géométrie).
- Le point le plus proche et les intersections sont calculés segment par segment (O(nombre de segments)).
- La méthode d'ajout de points à la souris trie les points ; deux points d'abscisses presque égales sont refusés.


