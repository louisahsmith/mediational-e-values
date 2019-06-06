# Computing code for *Mediational E-values: Approximate Sensitivity Analysis for Unmeasured Mediator-Outcome Confounding*
## To be published in *Epidemiology*

The results can be replicated using the R scripts in this repository. The file `functions.R` contains the functions used to generate the data as well as produce the figures and tables. The file `make.R` applies these functions in the proper sequence and populates the `eAppendix.Rmd` file, producing `eAppendix.pdf` (the tables and figures of current document). Only the `make.R` script needs to be run to replicate the results.

The size of the simulation can be controlled by the object `N`. Note that large `N` (e.g., greater than the 5,000,000 to produce the results in this paper) requires several hours of computing time and may overload a personal computer.

In the code, the variable names correspond to the following probabilities:

| name | probability                       |
|-----|---------------------------------|
| `a` | Pr(U = 1)                       |
| `b` | Pr(M = 1 \| A = 1, U = 1)        |
| `c` | Pr(M = 1 \| A = 1, U = 0)        |
| `d` | Pr(M = 1 \| A = 0, U = 1)        |
| `e` | Pr(M = 1 \| A = 0, U = 0)        |
| `f` | Pr(Y = 1 \| A = 1, U = 1, M = 1) |
| `g` | Pr(Y = 1 \| A = 1, U = 1, M = 0) |
| `h` | Pr(Y = 1 \| A = 1, U = 0, M = 1) |
| `i` | Pr(Y = 1 \| A = 1, U = 0, M = 0) |
| `j` | Pr(Y = 1 \| A = 0, U = 1, M = 1) |
| `k` | Pr(Y = 1 \| A = 0, U = 1, M = 0) |
| `l` | Pr(Y = 1 \| A = 0, U = 0, M = 1) |
| `m` | Pr(Y = 1 \| A = 0, U = 0, M = 0) |
| `n` | Pr(A = 1)                       |

For simplicity, we also computed the probabilities of each of A, M, U, and Y taking on values of 0. Since we assumed that each is binary, Pr(U = 0) = 1 - Pr(U = 1), Pr(M = 0 | A = 1, U = 1) = 1 - Pr(M = 1 | A = 1, U = 1), etc. In the code, those values are referred to as `a_comp`, `b_comp`, etc., for each of the above letters.

The alternative RRs are labeled as follows:

| name |       parameter      |
|------|----------------------|
| `RRau`                                                                                                                                                                                                                                                                                                                                                                                                             | RR<sub>AU \|M</sub> |
| `RRuy`                                                                                                                                                                                                                                                                                                                                                                                                             | RR<sub>UY \|(A = 1, M)</sub> |
| `RRam`                                                                                                                                                                                                                                                                                                                                                                                                             |  RR<sub>AM</sub>   |
| `RRam_no0`                                                                                                                                                                                                                                                                                                                                                                                                         |  RR<sub>AM=1</sub>  |
| `RRam_strat`                                                                                                                                                                                                                                                                                                                                                                                                       |   RR<sub>AM\|U</sub> |
| `RRam_strat_no0`                                                                                                                                                                                                                                                                                                                                                                                                   | RR<sub>AM=1\|U</sub> |
| `RRum`                                                                                                                                                                                                                                                                                                                                                                                                             | RR<sub>UM</sub> |
| `RRum_no0`                                                                                                                                                                                                                                                                                                                                                                                                         | RR<sub>UM=1</sub> |
| `RRum_strat`                                                                                                                                                                                                                                                                                                                                                                                                       | RR<sub>UM\|A</sub> |
| `RRum_strat_no0`                                                                                                                                                                                                                                                                                                                                                                                                   |  RR<sub>UM=1\|A</sub> |


After running the `make.R` script, the dataset containing all of the simulations can be loaded into the R environment with the command `drake::loadd(full_dat)`.
