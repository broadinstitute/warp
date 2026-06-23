### Glimpse2MergeBatches AF and INFO score recalculation
**By: Christopher Kachulis**

Glimpse outputs three values in the info field:

1) AF.  This is easily calculated for the full cohort based on the individual batch values as a weighted mean of the batch allele frequencies, $$AF_{cohort}=\frac{\sum AF_i N_i}{\sum N_i}$$

2) RAF.  This is the reference panel allele frequency.  Assuming the same reference panel was used for all batches, these are all the same, so just take the first value

3) INFO.  This is the "IMPUTE style info score".  Based on https://static-content.springer.com/esm/art%3A10.1038%2Fnrg2796/MediaObjects/41576_2010_BFnrg2796_MOESM3_ESM.pdf, this is calculated as $$1-\frac{\sum f_j - e^2_j}{2N\times AF(1-AF)}$$, or 1 if AF=0,1, with $j$ running across the N samples, $f_j = p_{j1} + 4 p_{j2}$, $e_j = p_{j1} + 2 p_{j2}$, where $p_{j1}$ is the imputed posterior of a het for sample $j$, and $p_{j2}$ is the imputed posterior of a hom var for sample $j$.  Note that the terms in the denominator are all easily calculated for a cohort based on their values for the constituent batches;  AF as described above, and N just as the sum over the batches.  The numerator we can define for batch $i$ as $C_i = \sum (f_{ij} -e_{ij}^2)$ and note that, for a whole cohort, we simply have $C_{cohort} = \sum C_i$.  We then note that, for a batch i, we have $$I_i =1 - \frac{C_i}{2N_i\times AF_i(1-AF_i)}$$, so we can solve for $C_i$ as $$C_i = (1-I_i)*2N_i\times AF_i(1-AF_i)$$.  We can then calculate the cohort INFO score as $$I_{cohort}=1-\frac{\sum C_i}{2N_{cohort} AF_{cohort}(1-AF_{cohort})}$$ which becomes, $$I_{cohort}=1-\frac{\sum (1-I_i)*2N_i\times AF_i(1-AF_i)}{2\sum N_i \times \frac{\sum AF_i N_i}{\sum N_i} (1-\frac{\sum AF_i N_i}{\sum N_i})}$$
