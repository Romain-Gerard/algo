import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

# ========== Factorielle et coefficient binomial ==========

def facto_rec(x):
    """
    Calcule la factorielle de x de manière récursive.
    """
    if x == 0:
        return 1
    return x * facto_rec(x-1)


def coeff_bino(n, k):
    """
    Calcule le coefficient binomial C(n, k).
    """
    return facto_rec(n) // (facto_rec(k)*facto_rec(n-k))


def coeff_bino_rec(n, k):
    """
    Calcule le coefficient binomial C(n, k) de manière récursive.
    """
    if k > n:
        return 0
    if k == 0 or k == n:
        return 1
    
    return coeff_bino_rec(n-1, k-1) + coeff_bino_rec(n-1, k)


def coeff_bino_dyn(n, k, memo={}):
    """
    Calcule le coefficient binomial C(n, k) de manière récursive et dynamique.
    """
    if k > n:
        return 0
    if k == 0 or k == n:
        return 1
    
    if (n, k) not in memo:
        memo[(n, k)] = coeff_bino_dyn(n-1, k-1) + coeff_bino_dyn(n-1, k)
    return memo[(n, k)]


# ========== Permutations et combinaisons ==========

def permute_rec(lst):
    """
    Retourne la liste de toutes les permutations de lst, de manière récursive.
    """
    if len(lst) == 0:
        return [[]]

    result = []
    for i in range(len(lst)):
        current = lst[i]
        remaining = lst[:i] + lst[i+1:]
        for e in permute_rec(remaining):
            result.append([current] + e)
    return result

def combine_dyn(lst, memo=None):
    """
    Génère toutes les combinaisons possibles de toutes les tailles à partir d'une liste.
    """
    if memo is None:
        memo = {}

    # Si la liste est vide, il n'y a qu'une seule combinaison : l'ensemble vide.
    if len(lst) == 0:
        return [[]]

    # Si le résultat est déjà calculé, on le retourne depuis la mémoire.
    lst_tuple = tuple(lst)
    if lst_tuple in memo:
        return memo[lst_tuple]

    # Récurrence :
    # On exclut le premier élément et génère les combinaisons du reste.
    combinaisons_sans_premier = combine_dyn(lst[1:], memo)
    # On inclut le premier élément dans chaque combinaison du reste.
    premier = lst[0]
    combinaisons_avec_premier = [[premier] + combinaison for combinaison in combinaisons_sans_premier]
    # On combine les deux ensembles de combinaisons.
    result = combinaisons_sans_premier + combinaisons_avec_premier

    memo[lst_tuple] = result

    return result


# ========== Nombre de Stirling de seconde espèce ==========

def calcule_stirling_rec(n, k):
    """
    Calcule le nombre de Stirling de seconde espèce S(n, k)
    selon une définition récursive.
    
    S(n, k) = S(n-1, k-1) + k * S(n-1, k)
    avec les cas de base :
       - S(0, 0) = 1
       - S(n, 0) = 0 pour n > 0
       - S(0, k) = 0 pour k > 0
    """
    if n == 0 and k == 0:
        return 1
    if (n == 0 and k > 0) or (n > 0 and k == 0):
        return 0

    return calcule_stirling_rec(n-1, k-1) + k*calcule_stirling_rec(n-1, k)


def calcule_stirling_dyn(n, k, memo={}):
    """
    Calcule le nombre de Stirling de seconde espèce (version dynamique).
    """
    if n == 0 and k == 0:
        return 1
    if (n == 0 and k > 0) or (n > 0 and k == 0):
        return 0

    if (n, k) not in memo:
        memo[(n, k)] = calcule_stirling_dyn(n-1, k-1) + k*calcule_stirling_dyn(n-1, k)
    return memo[(n, k)]


# ========== Nombre de Bell ==========

def calcule_bell_iter(n):
    res = 0
    for k in range(1, n+1):
        res += calcule_stirling_rec(n, k)
    return res


def calcule_bell_rec(n):
    if n == 0:
        return 1
    
    res = 0
    for m in range(n):
        res += coeff_bino_rec(n-1, m) * calcule_bell_rec(m)
    return res


def calcule_bell_dyn(n, memo={}):
    # Version dynamique.
    if n == 0:
        return 1
    
    if n in memo:
        return memo[n]
    
    res = 0
    for m in range(n):
        res += coeff_bino_dyn(n-1, m) * calcule_bell_dyn(m)

    memo[n] = res
    return res


# ========== Partitions de type A ==========

def partitionne_stirling_rec(lst, k):
    """
    Retourne la liste de toutes les partitions de lst en k sous-listes (blocs) en
    s"inspirant de la définition récursive du nombre de Stirling de seconde espèce:
    
                        S(n, k) = S(n-1, k-1) + k * S(n-1, k)
    """
    n = len(lst)

    # -------------
    # Cas de base :
    # -------------

    # Si k = 0 et n = 0, alors on a une seule partition vide.
    if k == 0 and n == 0:
        return [[]]
    
    # Si k = 0 ou n = 0, mais pas les deux en même temps,
    # alors il n"existe aucune partition.
    if (k == 0 or n == 0) and (k != n):
        return []
    
    # Si k = n, chaque élément de lst doit être dans son propre bloc.
    if k == n:
        return [ [[x] for x in lst] ]
    
    # Si k = 1, il n"y a qu"un seul bloc contenant tous les éléments.
    if k == 1:
        return [[lst]]
    
    # Si k > n, il n"existe pas de partition car on ne peut
    # pas répartir n éléments en plus de n blocs non vides.
    if k > n:
        return []
    
    # ----------------------------------------------------
    # Récurrence :
    # 1) On met lst[0] dans son propre bloc, et on partitionne
    #    lst[1:] en (k-1) blocs.
    # 2) On met lst[0] dans chacun des blocs d"une partition de
    #    lst[1:] en k blocs.
    # ----------------------------------------------------

    # Initialiser la liste de toutes les partitions.
    lst_part = []
    
    # 1) Ajouter lst[0] dans un nouveau bloc (S(n-1, k-1))
    for part in partitionne_stirling_rec(lst[1:], k - 1):
        # On ajoute le bloc [lst[0]] devant
        new_part = [[lst[0]]] + part
        lst_part.append(new_part)
    
    # 2) Ajouter lst[0] dans chacun des blocs.
    for part in partitionne_stirling_rec(lst[1:], k):
        # Pour chaque partition obtenue, on crée k nouvelles partitions
        # où lst[0] est inséré dans un seul bloc par partition.
        for i in range(len(part)):
            new_part = []
            for j, bloc in enumerate(part):
                if j == i:
                    # On ajoute lst[0] dans le bloc d"indice i ...
                    new_part.append([lst[0]] + bloc)
                else:
                    # ... sinon on recopie le bloc.
                    new_part.append(bloc[:])
            lst_part.append(new_part)
    
    return lst_part


def partitionne_stirling_dyn(lst, k, memo=None):
    """
    Version dynamique de partitionne_stirling_rec qui utilise la mémorisation.
    """
    if memo is None:
        memo = {}
    
    n = len(lst)
    
    # Créer une clé unique
    key = (tuple(lst), k)
    
    # Si le résultat est déjà dans la mémoire, le retourner.
    if key in memo:
        return memo[key]

    # -------------
    # Cas de base :
    # -------------
    if k == 0 and n == 0:
        return [[]]
    
    if (k == 0 or n == 0) and (k != n):
        return []
    
    if k == n:
        return [ [[x] for x in lst] ]
    
    if k == 1:
        return [[lst]]
    
    if k > n:
        return []

    # Initialiser la liste des partitions
    lst_part = []
    
    # 1) Ajouter lst[0] dans un nouveau bloc
    for part in partitionne_stirling_dyn(lst[:-1], k - 1, memo):
        new_part = part + [[lst[-1]]]
        lst_part.append(new_part)
    
    # 2) Ajouter lst[0] dans chacun des blocs existants
    for part in partitionne_stirling_dyn(lst[:-1], k, memo):
        for i in range(len(part)):
            new_part = []
            for j, bloc in enumerate(part):
                if j == i:
                    new_part.append(bloc + [lst[-1]])
                else:
                    new_part.append(bloc[:])
            lst_part.append(new_part)
    
    # Mémoriser le résultat.
    memo[key] = lst_part
    
    return lst_part

def partitionne_bell(lst):
    """
    Retourne la liste de toutes les partitions des éléments de lst en s"inspirant
    de la définition du nombre de Bell comme une somme de nombres de Stirling de
    seconde espèce.
    """
    n = len(lst)
    all_part = []
    for k in range(n+1):
        all_part = all_part + partitionne_stirling_dyn(lst, k)
    return all_part


# ========== Statistiques sur les permutations de type A ==========

def liste_desc(perm):
    """
    Retourne les indices i tels que σ(i) > σ(i+1).
    """
    n = len(perm)
    lst = [i for i in range(n-1) if perm[i] > perm[i+1]]
    return lst


def liste_asc(perm):
    """
    Retourne les indices i tels que σ(i) < σ(i+1).
    """
    n = len(perm)
    lst = [i for i in range(n-1) if perm[i] < perm[i+1]]
    return lst


def liste_exce(perm):
    """
    Retourne les indices i tels que σ(i) > i.
    """
    n = len(perm)
    lst = [i for i in range(n-1) if perm[i] > i]
    return lst


def liste_inv(perm):
    """
    Retourne les paires (i, j) avec i < j et σ(i) > σ(j).
    """
    n = len(perm)
    lst = [(i, j) for i in range(n) for j in range(i+1, n) if perm[i] > perm[j]]
    return lst


# ========== Partitions de type B ==========

# Une partition de type B est une partition d'un ensemble ⟨n⟩ qui respecte les propriétés suivantes :
# 1. Pour tout i ≥ 1, les blocs π2i et π2i−1 sont opposés, c'est-à-dire que π2i = −π2i−1,
#    où −β représente l'ensemble des opposés des éléments de β (−β = {−a : a ∈ β}).
# 2. Il existe un bloc spécial appelé le bloc zéro (π0), qui est symétrique, c'est-à-dire que si a ∈ π0,
#    alors −a ∈ π0 également.
#
# Représentation d'Adler :
# - Dans le bloc zéro, tous les éléments négatifs sont supprimés.
# - Dans chaque bloc, les éléments sont triés par valeur absolue croissante.
# - Pour chaque paire de blocs opposés, seul le bloc avec le premier élément positif est conservé.

# Stratégie pour générer une partition de type B d'un ensemble de taille et divisée en k blocs :
# 1.1. Récupérer une partition de type A d'un ensemble de taille n divisée en k blocs (liste de k sous-listes).
# 1.2. Trier tous les éléments de chaque sous-liste par valeur croissante.
# 1.3. Trier les sous-listes par valeur croissante de leur premier élément.
# 1.4. Récupérer la liste des indices des éléments qui ne sont pas les premiers de leur bloc et pas dans le bloc zéro.
# 2. Générer toutes les combinaisons des indices de la liste obtenue à l'étape 1.4. de longueur 0 à la longueur de la liste.
# 3. Pour chaque combinaison d'indices, créer une nouvelle partition de type B en inversant le signe aux indices de la combinaison.


def print_partitions(parts, length=-1, name="partitions", scope=5):
    """
    Affiche les partitions de manière lisible.
    """
    if len(parts) > 2 * scope:
        for i in range(scope):
            print(parts[i])
        print("...")
        for i in range(-scope, 0):
            print(parts[i])
    else:
        for i in range(len(parts)):
            print(parts[i])
    print(f"Nombre de {name} pour n={length} : {len(parts)}")


def sort_partitions(partitions):
    """
    Trie les partitions de type B (ou A) selon la représentation d'Adler.
    """
    sorted_partitions = []
    for part in partitions:
        # Trier les blocs par valeur croissante de leur premier élément
        part.sort(key=lambda x: x[0])
        # Trier les éléments de chaque bloc par valeur absolue croissante
        for bloc in part:
            bloc.sort(key=lambda x: abs(x))
        sorted_partitions.append(part)
    return sorted_partitions


def get_signable_inds(part):
    """
    Récupère la liste des indices des éléments qui ne sont pas les premiers de leur bloc et pas dans le bloc zéro.
    """
    liste_signable_inds = []
    ind = 0
    for i in range(len(part)):
        for j in range(0, len(part[i])):
            if part[i][0] != 0 and part[i][j] != part[i][0]:
                liste_signable_inds.append(ind)
            ind += 1
    return liste_signable_inds


def convert_part_A_to_B(part_A, permu_signed_inds):
    """
    Convertit une partition de type A en une partition de type B en inversant le signe
    aux indices spécifiés par permu_signed_inds.
    """
    ind = 0
    part_B = []
    for i in range(len(part_A)):
        bloc = []
        for j in range(0, len(part_A[i])):
            if ind in permu_signed_inds:
                bloc.append(-part_A[i][j])
            else:
                bloc.append(part_A[i][j])
            ind += 1
        part_B.append(bloc)
    return part_B


def get_all_parts_B_from_A(parts_A, complete = False):
    """
    Récupère toutes les partitions de type B à partir d'une liste de partitions de type A.
    """
    all_part_B = []
    for part_A in parts_A:
        # Trier la partition de type A selon la représentation d'Adler
        #part_A = sort_partitions([part_A])[0]
        # Récupérer la liste des indices signables
        liste_signable_inds = get_signable_inds(part_A)
        # Générer toutes les combinaisons d'indices signables
        combinaisons = combine_dyn(liste_signable_inds)
        for comb in combinaisons:
            # Convertir la partition de type A en type B
            part_B = convert_part_A_to_B(part_A, comb)
            all_part_B.append(part_B)

    if complete:
        all_part_B = complete_parts_B(all_part_B)

    return all_part_B

def complete_parts_B(parts_B):
    """
    Complète les partitions de type B en ajoutant après charque bloc son bloc opposé.
    """
    all_part_B = []
    for part_B in parts_B:
        new_part = []
        for bloc in part_B:
            if 0 not in bloc:
                new_part.append([-x for x in bloc])
                new_part.append(bloc)
            else:
                new_part.append(bloc + [-x for x in bloc if x != 0])
        all_part_B.append(new_part)
    return all_part_B

def calcule_dowling(n):
    """
    Fonction pour calculer le nombre de partitions de type B selon la formule du nombre de Dowling
    """
    dn = 0
    for i in range(n+1):
        # coeff binomial C(n, i)
        c = coeff_bino_dyn(n, i)
        # somme sur k
        s = 0
        for k in range(n - i + 1):
            s += (2**(n - i - k)) * calcule_stirling_dyn(n - i, k)
        dn += c * s
    return dn


def calcule_dowling_no_zero_block(n):
    """
    Fonction pour calculer le nombre de partitions de type B sans bloc zéro selon la formule du nombre de Dowling
    """
    wn = 0
    for k in range(0, n+1):
        wn += (2**(n - k)) * calcule_stirling_dyn(n, k)
    return wn


# ========== Type B separated set partitions =========

def is_separated(part):
    """
    Vérifie si une partition de type B est séparée.
    """
    # Vérifier les blocs
    for bloc in part:
        # Vérifier si les éléments sont consécutifs
        for i in range(len(bloc) - 1):
            if abs(bloc[i]) + 1 == abs(bloc[i + 1]) and bloc[i] * bloc[i + 1] > 0:
                return False
    return True

def is_strongly_separated(part):
    """
    Vérifie si une partition de type B est fortement séparée.
    """
    # Vérifier les blocs
    for bloc in part:
        # Vérifier si les éléments sont consécutifs
        for i in range(len(bloc) - 1):
            if abs(bloc[i]) + 1 == abs(bloc[i + 1]):
                return False
    return True


# ========== Type B merging-free partitions =========

def is_merge_free(part):
    """
    Vérifie si une partition de type B est merge-free.
    """
    for i in range(1, len(part)):
        # Vérifier si le maximum du bloc i-1 est inférieur au minimum du bloc i en valeur absolue
        if max(abs(np.array(part[i-1]))) < min(abs(np.array(part[i]))):
            return False
    return True

def is_normal_merge_free(part):
    """
    Vérifie si une partition de type B est merge-free.
    """
    for i in range(1, len(part)):
        # Vérifier si le maximum du bloc i-1 est inférieur au minimum du bloc i.
        if max(part[i-1]) < min(part[i]):
            return False
    return True


# ========== Statistiques sur les partitions de type B ==========

def compte_inversions(part):
    """
    Compte le nombre d'inversions dans une partition.
    """
    n = len(part)
    # On calcule le minimum de chaque bloc en valeur absolue
    lst_min = np.array([abs(np.array(bloc)).min() for bloc in part])
    # On parcourt tous les blocs
    inv_count = 0
    for i in range(n-1):
        # On parcourt tous les éléments du bloc i
        for j in range(len(part[i])):
            # On vérifie si l'élément est supérieur à au moins un des minimums des blocs suivants
            test = part[i][j] > lst_min[i+1:]
            inv_count += sum(test)
    return inv_count


# ========== Type B inversion free set partitions =========

def is_inversion_free(part):
    """
    Vérifie si une partition de type B est inversion-free
    """
    n = len(part)
    # On calcule le minimum de chaque bloc en valeur absolue
    lst_min = np.array([abs(np.array(bloc)).min() for bloc in part])
    # On parcourt tous les blocs
    for i in range(n-1):
        # On parcourt tous les éléments du bloc i
        for j in range(len(part[i])):
            # On vérifie si l'élément est supérieur à au moins un des minimums des blocs suivants
            test = part[i][j] > lst_min[i+1:]
            if True in test:
                return False
    return True


# ========== Type B non-nesting partitions =========

def plot_partition(partition, figsize=(10, 2), arc_height_scale=1.0):
    """
    Affiche une partition d'ensemble sous forme d'arcs, en ordonnant d'abord
    toutes les valeurs positives croissantes, puis toutes les valeurs négatives
    par ordre d'absolu croissant (ex. 1, 2, 3, -1, -2, -3).

    Args:
        partition (list of list of int]): blocs contenant des entiers positifs et négatifs.
        figsize (tuple): taille de la figure matplotlib (largeur, hauteur).
        arc_height_scale (float): facteur d'échelle pour la hauteur des arcs.
    """
    # 1) Calcule l'ordre des éléments : (groupe, clé)
    #    groupe = 0 pour x >= 0, 1 pour x < 0 ; clé = abs(x)
    elements = sorted(
        {x for bloc in partition for x in bloc},
        key=lambda x: (x < 0, abs(x))
    )

    # 2) Mappe chaque valeur à une position 1,2,…,n
    pos = {x: i + 1 for i, x in enumerate(elements)}
    n = len(elements)

    fig, ax = plt.subplots(figsize=figsize)

    # 3) Trace les points et leurs labels
    for x in elements:
        ax.plot(pos[x], 0, 'o', color='k')
        ax.text(pos[x], -0.1, str(x), ha='center', va='top')

    # 4) Pour chaque bloc, dessine les arcs entre ses voisins dans ce nouvel ordre
    for bloc in partition:
        # tri du bloc selon la même clé
        tri_bloc = sorted(bloc, key=lambda x: (x < 0, abs(x)))
        for i in range(len(tri_bloc) - 1):
            v1, v2 = tri_bloc[i], tri_bloc[i + 1]
            x1, x2 = pos[v1], pos[v2]
            mid = 0.5 * (x1 + x2)
            width = x2 - x1
            height = width * arc_height_scale
            arc = Arc(
                (mid, 0),
                width=width, height=height,
                angle=0, theta1=0, theta2=180,
                edgecolor='black'
            )
            ax.add_patch(arc)

    # 5) Ajuste les limites et masque les axes
    ax.set_xlim(0, n + 1)
    ax.set_ylim(-1, n * arc_height_scale / 2 + 1)
    ax.axis('off')
    plt.tight_layout()
    plt.show()


# ========== Stirling Permutations ==========

def get_tree_from_stir_perm(perm):

    visited = []
    graph = []

    for e in perm:
        if len(visited) == 0:
            graph.append((0, e))
            visited.append(e)
        else:
            if e not in visited:
                graph.append((visited[-1], e))
                visited.append(e)
            else:
                visited = visited[:-1]

    return graph

def count_leaves_on_stir_perm(perm):
    count = 0
    for i in range(len(perm) - 1):
        if perm[i] == perm[i+1]:
            count += 1
    return count

def generate_stirling_permutations_rec(n):
    """
    Version récursive.
    """
    # Cas de base
    if n == 1:
        return [[1, 1]]
    
    all_perm = []
    for perm in generate_stirling_permutations_rec(n - 1):
        # On insère nn à chaque position possible
        for i in range(len(perm) + 1):
            new_perm = perm[:i] + 2*[n] + perm[i:]
            all_perm.append(new_perm)

    return all_perm


# ========== Flattened Stirling Permutations ==========

def is_run_sorted(perm):
    """
    Vérifie si une permutation est run-sorted.
    """
    min = perm[0]
    for i in range(len(perm) - 1):
        if perm[i] > perm[i+1]:
            if perm[i+1] >= min:
                min = perm[i+1]
            else:
                return False
    return True


def get_flattened_stirling_permutations(part_B):
    """
    On traduit une partition de type B en une permutation de Stirling.
    """
    if len(part_B) == 0:
        return []

    # If a in absolute value > min(part_B[i+1]), place the barred elements before.
    for i in range(len(part_B) - 1):
        block_inf = [e for e in part_B[i] if abs(e) < part_B[i+1][0]]
        block_sup_neg = [e for e in part_B[i] if (abs(e) > part_B[i+1][0]) * (e < 0)]
        block_sup_pos = [e for e in part_B[i] if (abs(e) > part_B[i+1][0]) * (e > 0)]
        part_B[i] = block_inf + block_sup_neg + block_sup_pos

    # Copie profonde
    part = [block[:] for block in part_B]

    # 1)
    # Remplacer 1 par 1 1
    part[0] = [1] + part[0]

    # 2)
    # If a is in the i-th block with a < m_{i+1}, i < k, or in the last block, and it is barred, then insert it between the 1’s.
    for i in range(len(part_B) - 1):
        for j in range(len(part_B[i])):
            if part_B[i][j] < 0 and abs(part_B[i][j]) < part_B[i+1][0]:
                index = max(i for i, x in enumerate(part[0]) if x == 1)
                part[i].remove(part_B[i][j])
                part[0].insert(index, part_B[i][j])
    # In the last block, insert the barred elements between the 1s.
    for j in range(len(part_B[-1])):
        if part_B[-1][j] < 0:
            index = max(i for i, x in enumerate(part[0]) if x == 1)
            part[-1].remove(part_B[-1][j])
            part[0].insert(index, part_B[-1][j])

    # 3)
    # For each i >= 1, insert m_{i+1} in the i-th block after the barred elements.
    for i in range(len(part) - 1):
        min_next = part[i+1][0]
        if any(x < 0 and abs(x) > min_next for x in part[i]):
            ind_last_barred = max(i for i, x in enumerate(part[i]) if x < 0)
            part[i].insert(ind_last_barred + 1, min_next)
        else:
            ind = max(i for i, x in enumerate(part[i]) if abs(x) < min_next)
            part[i].insert(ind + 1, min_next)

    # 4)
    # Replace each element a different from the minimum with aa.
    # 5)
    # Remove the bars and the separators.
    perm = []
    for i in range(len(part)):
        min_current = part[i][0]
        if i < len(part) - 1:
            min_next = part[i+1][0]
        else:
            min_next = float('inf')
        for j in range(len(part[i])):
            if part[i][j] != min_current and part[i][j] != min_next:
                perm.append(part[i][j])
                perm.append(part[i][j])
            else:
                perm.append(part[i][j])
    perm = [abs(x) for x in perm]

    return perm
