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
            new_part.append([-x for x in bloc])
            new_part.append(bloc)
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

