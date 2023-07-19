[**Status:** Archive (code is provided as-is, no updates expected)]: #
# CNMAC2023

This repo contains the code for implementing example 3.1 of paper [*Invariance Principle and the Asymptotic Behavior of T-S Fuzzy Systems*](CNMAC_2023_Fuzzy.pdf).

## Contents

This repo contains
| | |
|:---|:---:|
|README.md|This file|
|example3_1.m|Example 3.1 script|
|CNMAC_2023_Fuzzy.pdf|Preprint .pdf version of the paper|
|LICENSE|GNU General Public License v3.0|

## Requirements

- Matlab R2021b
- [Sedumi](https://sedumi.ie.lehigh.edu/)
- [Yalmip](https://yalmip.github.io/)

## How to use it

Open the `example3_1.m` file in matlab and run the script.
The three figures below should open.

![Figure 1](/figs/fig1.bmp)

Figure 1 shows the graph of $V(x)$ inside set $Z$.  The minimum value of $V$ in the boundary of $Z$ is $b$.

![Figure 2](/figs/fig2.bmp)

Figure 2 shows the contour plot of $V$ superposed by the boundary of set $\mathcal{D}$, in black. Using the level sets, we can find $l$, the maximum value of $V$ inside $\mathcal{D}$.

![Figure 3](/figs/fig3.bmp)

Figure 3 shows the level sets $\Omega_L,\Omega_l$ and the set $\mathcal{D}$.

## Citation

If you find this paper useful, please cite it in your publication using the following bibtex entry:

```
@inproceedings{araujoInvariancePrincipleAsymptotic2023,
  title = {Invariance Principle and the Asymptotic Behavior of T-S Fuzzy Systems},
  booktitle = {Proceeding Series of the Brazilian Society of Computational and Applied Mathematics},
  author = {Araujo, Rayza and Alberto, Luis Fernando Costa and Valentino, Michele Cristina},
  year = {2023},
  address = {{Bonito, MS, Brazil.}}
}
```