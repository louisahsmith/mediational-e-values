# Computing code for *Mediational E-values: Approximate Sensitivity Analysis for Unmeasured Mediator-Outcome Confounding*, to be published in *Epidemiology*

The results can be replicated using the R scripts in the supplementary material. The file `functions.R` contains the functions used to generate the data as well as produce the figures and tables. The file `make.R` applies these functions in the proper sequence and populates the `eAppendix.Rmd` file, producing `eAppendix.pdf` (the tables and figures of current document). Only the `make.R` script needs to be run to replicate the results.

The size of the simulation can be controlled by the object `N`. Note that large `N` (e.g., greater than the 5,000,000 to produce the results in this paper) requires several hours of computing time and may overload a personal computer.

In the code, the variable names correspond to the following probabilities:

| | |
|`a`| = Pr(A = 1)|
|`b`| = Pr(U = 1)|
|`c`| = |
|`d`| = |
|`e`| = |
|`f`| = |
|`g`| = |
|`h`| = |
|`i`| = |
|`j`| = |
|`k`| = |
|`l`| = |
|`m`| = |
|`n`| = |

(Strictly speaking, the bias we are interested in does not depend on $\pr(Y = 1 \mid A = 0, M = m, U= u)$ for any $m$ or $u$, but we included those values in order to calculate the true and observed NDEs and not just their ratio. The marginal probability of the exposure, $\Na$, is only useful for the calculation of one of the possible RR alternatives and also does not directly factor into the bias or the effect sizes.)

For simplicity, we also computed the probabilities of each of $A$, $M$, $U$, and $Y$ taking on values of 0. Since we assumed that each is binary, $\pr(U = 0) = 1 - \Au$, $\pr(M = 0 \mid A = 1, U = 1) = 1 - \Bm$, etc. In the code, those values are referred to as `a_comp`, `b_comp`, etc., for each of the above letters.

The alternative RRs are labeled as follows:

\texttt{RRau} & = \RRau \\
\texttt{RRuy} & = \RRuy \\
\texttt{RRam} & = \text{RR}_{AM} \\
\texttt{RRam\_no0} & = \text{RR}_{AM=1} \\
\texttt{RRam\_strat} & = \text{RR}_{AM|U} \\
\texttt{RRam\_strat\_no0} & = \text{RR}_{AM=1|U} \\
\texttt{RRum} & = \text{RR}_{UM} \\
\texttt{RRum\_no0} & = \text{RR}_{UM=1} \\
\texttt{RRum\_strat} & = \text{RR}_{UM|A} \\
\texttt{RRum\_strat\_no0} & = \text{RR}_{UM=1|A}


After running the `make.R` script, the dataset containing all of the simulations can be loaded into the R environment with the command `drake::loadd(full_dat)`.
