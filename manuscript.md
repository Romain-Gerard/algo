## Type B strongly separated set partitions
A strongly separated partition is a partition in which all integers within the same block are non-consecutive, regardless of their sign.

We give a bijection between type B set partitions without zero block over $\langle n \rangle$ and type B strongly separated set partitions whithout zero block over \langle n+1 \rangle.

Let $SS^B_n$ denote the set of type B strongly separated partitions without zero block over \<n>.

We define the map $g\ : \Pi^0_n \mapsto SS^B_{n+1}$ by $g(\pi) = \pi'$, where $\pi'$ is obtained from $\pi$ as follows.

Let $\pi = \pi_1 | \pi_2 | \cdots | \pi_k \in \Pi^0_n$. We apply the following procedure to $\pi$.

1. Add a zero block which only contains $0$.
2. For $a$ from $2$ to $n$, if $a$ (or $\overline{a}$) is a succession in block $\pi_i$, $i>1$, then move $a$ (or $\overline{a}$) to the zero block.
3. Increase every integer in absolute value by $1$.

**Example 1.** Let $\pi = 1\ 4|2\ 6\ \overline{7}|3\ 8\ 9\ 10 \ \in \Pi^0_{10}$. Then, we add a zero block. Since $\overline{7}$ and $9$ are successions, we move theme to the zero block. Consequently, $10$ is no longer a succession. Then we increase every integer in absolute value by 1. Therefore, $\pi' = 1\ \overline{8}\ 10 | 2\ 5 | 3\ 7 | 4\ 9\ 10\ \in SS^B_{11}$.

Let us define $\pi = g^{-1}(\pi')$ as follows. Let $\pi' = \pi_1' | \pi_2' | \cdots | \pi_k' \in SS^B_{n+1}$.

1. Decrease every integer in absolute value by 1.
2. For every element $a$ (or $\overline{a}$) $\ne 0$ in the first block, move $a$ (or $\overline{a}$) to the block containing $a-1$ (or $\overline{a-1}$).
3. Remove the first block, which has only the element $0$ left.

**Example 2.** Let $\pi' = 1\ \overline{8}\ 10 | 2\ 5 | 3\ 7 | 4\ 9\ 10\ \in SS^B_{11}$. Then, we decrease every integer in absolute value by 1. Now, we have $\pi_1 = 0\ \overline{7}\ 9$. Then, we move $\overline{7}$ to the third block and $9$ to the fourth block. Finally, we remove $\pi_1$. Therefore, $\pi = 1\ 4|2\ 6\ \overline{7}|3\ 8\ 9\ 10 \ \in \Pi^0_{10}$.

We note that

**Lemma.** *For any partition* $\pi \in SS^B_n$ *, we have* $(g \circ g^{-1})(\pi) = \pi$.

**Proposition.** *The map g is a bijection.*
