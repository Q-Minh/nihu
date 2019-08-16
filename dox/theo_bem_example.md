An introductory BEM example {#theo_bem_example}
===========================

\page theo_bem_example An introductory BEM example

[TOC]

Problem definition {#bem_example_def}
==================

Let us consider the Helmholtz equation with a Neumann boundary condition:

\f$ 
\displaystyle
\left(\nabla^2 + k^2 \right) p({\bf x}) = 0 \qquad {\bf x} \in V \\
\displaystyle
\partial p({\bf x}) / \partial n = q({\bf x}) = \bar{q}({\bf x}) \qquad {\bf x} \in S
\f$ 

where  \f$ V\subset\mathbb{R}^d \f$ is a bounded domain with a smooth boundary  \f$ S \f$ ,  \f$ p \f$ denotes the acoustic pressure,  \f$ q \f$ denotes its normal derivative (velocity),  \f$ \bar{q} \f$ denotes the prescribed normal velocity on the boundary  \f$ S \f$ , and  \f$ k \f$ denotes the acoustic wave number.
The pressure field  \f$ p({\bf x}) \f$ is sought on the boundary and in internal points of the domain  \f$ V \f$ .


Fundamental solution and integral operators {#bem_example_op}
===========================================

Solution of the boundary value problem is formulated using the fundamental solution  \f$ G({\bf x}_0, {\bf x}) \f$ of the Helmholtz equation defined by

\f$
\displaystyle
\left(\nabla^2 + k^2 \right) G({\bf x}_0, {\bf x}) = -\delta({\bf x}-{\bf x}_0), \qquad {\bf x} \in \mathbb{R}^d \setminus \left\{ {\bf x}_0 \right\}
\f$

where  \f$ \delta \f$ denotes the Dirac delta function.
The fundamental solution and its normal derivative are used as kernel functions to form the following single and double layer potential integral operators, respectively

\f$ 
\displaystyle
\left(\mathcal{L}u\right)_S({\bf x}) = \int_S G({\bf x},{\bf y}) u({\bf y}) \mathrm{d} S_{{\bf y}} \\
\displaystyle
\left(\mathcal{M}u\right)_S({\bf x}) = \int_S \frac{\partial G({\bf x}, {\bf y})}{\partial n_{{\bf y}}}({\bf x},{\bf y}) u({\bf y}) \mathrm{d} S_{{\bf y}}
\f$ 

The introduced operators are applied to pose an integral equation problem equivalent to the Helmholtz problem:

\f$ 
\displaystyle
\left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}) = \frac{1}{2} p({\bf x}) \qquad {\bf x} \in S \\
\displaystyle
\left(\mathcal{M}p\right)_S({\bf x}) - \left(\mathcal{L}q\right)_S({\bf x}) = p({\bf x}) \qquad {\bf x} \in V \setminus S
\f$ 

The first, boundary integral equation needs to be solved, using the prescribed velocity field  \f$ \bar{q}({\bf x}) \f$ 
as excitation, to obtain the acoustic pressure  \f$ p({\bf x}) \f$ on the boundary  \f$ S \f$ .
Consecutively, the acoustic pressure field  \f$ p({\bf x}) \f$ in  \f$ V \f$ can be
determined by evaluating the second equation.

Weighted residuals {#bem_example_res}
==================

The solution of the boundary integral equation is based on the weighted
residual method: The boundary integral equation is assumed to hold in a weak sense:

\f$ 
\displaystyle
\left< v, \left(\mathcal{M}p\right)_S - \left(\mathcal{L}q\right)_S - \frac{1}{2}p \right>_S = 0
\f$ 

where  \f$ \left<\cdot,\cdot\right>_S \f$ denotes the inner product defined as

\f$ 
\displaystyle
\left<v,f\right>_S = \int_{S} v({\bf x}) f({\bf x}) \mathrm{d} S_{{\bf x}}
\f$ 

and \f$ v({\bf x}) \f$ denotes an arbitrary test function chosen from the test function space  \f$ T \subset S \to \mathbb{R} \f$ .


Discretisation {#bem_example_disc}
==============

Numerical solution of the weak form is found by discretising the test space  \f$ T \f$ and the domain spaces of the boundary integral operators  \f$ \mathcal{L} \f$ and  \f$ \mathcal{M} \f$ .
Discretisation is performed by introducing their finite sets of basis functions  \f$ t_i \f$ ,  \f$ d^\mathcal{L}_j \f$ and  \f$ d^\mathcal{M}_j \f$ :

\f$ 
\displaystyle
v \in T \to v({\bf x}) = \sum_{i} v_i t_i({\bf x}) \\
\displaystyle
p \in D_{\mathcal{M}} \to p({\bf y}) = \sum_{j} p_j d^{\mathcal{M}}_j({\bf y}) \\
\displaystyle
q \in D_{\mathcal{L}} \to q({\bf y}) = \sum_{l} q_l d^{\mathcal{L}}_l({\bf y})
\f$ 

The discretised weak form of the boundary integral equation simplifies to a linear system of equations

\f$ 
\displaystyle
\sum_j \left< t_i, \left(\mathcal{M}d^{\mathcal{M}}_j\right)_S \right>_S p_j
-
\sum_l \left< t_i, \left(\mathcal{L}d^{\mathcal{L}}_l\right)_S \right>_S q_l
=
\frac{1}{2}
\left<t_i, d^{\mathcal{M}}_i\right>_S p_i
\f$ 

or, with matrix-vector notations,

\f$ 
\displaystyle
{\bf M} {\bf p} - {\bf L} {\bf q} = \frac{1}{2} {\bf B} {\bf p}
\f$ 

where  \f$ {\bf M} \f$ is a square dense matrix,  \f$ {\bf L} \f$ is a not necessarily square dense matrix and  \f$ {\bf B} \f$ is a sparse square matrix.
The Neumann problem is then solved by rearranging and inverting the linear equation

\f$ 
\displaystyle
{\bf p} = \left( {\bf M} - \frac{1}{2} {\bf B} \right)^{-1} {\bf L} \bar{{\bf q}}
\f$ 

BEM formalisms {#bem_example_form}
==============

Different BEM formalisms can be defined by classifying the selection of the function spaces.
- Galerkin BEM methods refer to the selection case, where the test function space coincides with the domain space of the solution:  \f$ t_i = d_i^{\mathcal{M}} \f$ .
In this selection case, the inner products are evaluated by double surface integrals, and, applying fem terminology, matrix  \f$ \bf B \f$ simplifies to a ,,boundary mass matrix''.
- Collocational BEM methods denote the case, where the test function space is generated by Dirac delta functions located at nodal locations  \f$ {\bf x}_i \f$ of the solution:  \f$ t_i({\bf x}) = \delta_{{\bf x}_i} =  \delta\left({{\bf x}-{\bf x}_i}\right) \f$ .
In the collocational case, the inner products boils down to single surface integrals, and matrix  \f$ \bf B \f$ typically simplifies to the unity matrix.

After the solution has been performed by means of direct or iterative techniques, the radiated field  \f$ p({\bf x}_i), {\bf x}_i \in V \f$ can be computed by evaluating the discretised version of the radiation integral

\f$ 
\displaystyle
\sum_j \left< \delta_{{\bf x}_i}, \left(\mathcal{M}d^{\mathcal{M}}_j\right)_S \right>_F p_j
-
\sum_k \left< \delta_{{\bf x}_i}, \left(\mathcal{L}d^{\mathcal{L}}_l\right)_S \right>_F q_l
=
p({\bf x}_i)
\f$ 


