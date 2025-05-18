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

