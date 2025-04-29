## Type B strongly separated set partitions
A strongly separated partition is a partition in which all integers within the same block are non-consecutive, regardless of their sign.

We give a bijection between type B set partitions without zero block over \<n> and type B strongly separated set partitions whithout zero block over <n+1>.

Let $SS^B_n$ denote the set of type B strongly separated partitions without zero block over <n>.

We define the map $g\ : \Pi^0_n \mapsto SS^B_{n+1}$ by $g(\pi) = \pi'$, where $\pi'$ is obtained from $\pi$ as follows.

Let $\pi = \pi_1 | \pi_2 | \cdots | \pi_k$. We apply the following procedure to $\pi$.

1. Add a zero block which only contains $0$.
2. For $a$ from $2$ to $n$, if $a$ (or $\overline{a}$) is a succession in block $\pi_i$, $i>1$, then move $a$ (or $\overline{a}$) to the zero block.
3. Increase every integer in absolute value by $1$.

**Exemple.** Let $\pi = 1\ 4|2\ 6\ \overline{7}$

