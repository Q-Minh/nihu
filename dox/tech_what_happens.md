What happens when NiHu evaluates a weighted residual? {#tech_what_happens}
=====================================================

\page tech_what_happens What happens when NiHu evaluates a weighted residual?

[TOC]

Introduction {#tech_what_happens_intro}
============

This page explains what happens when NiHu evaluates a weighted double integral

\f$
\displaystyle W_{ij} = \left< w_i, \left(\mathcal{K} d_j\right)_S \right>_F
\f$

or, as coded in a NiHu equivalent form
~~~~~~~~~~~
W << w * K[d];
~~~~~~~~~~~
where
- `w` and `d` are function spaces created either using function NiHu::create_function_space or using a view generating function like NiHu::constant_view, NiHu::isoparametric_view or NiHu::dirac.
- `K` is an integral operator created using function NiHu::create_integral_operator from a kernel function.
- and `W` is a matrix object

Operators and proxy objects {#tech_what_happens_operators}
===========================

When evaluating the exemplified expression, the integral transform (`[]` "index") operator NiHu::integral_operator::operator[] is evaluated first with the test function space as argument. This function simply returns a temporary NiHu::integral_transform instance. An NiHu::integral_transform is a proxy class that simply stores references to the NiHu::integral_operator `K` and the function_space `d`. The proxy class provides new functionalities: For example an NiHu::integral_transform can be tested by a test NiHu::function_space, as shown below.

The next step is the evaluation of the testing operator (`*` "multiplication") with a left hand side NiHu::function_space `w` and a right hand side NiHu::integral_transform (temporary). This testing operation returns a new proxy object of type NiHu::weighted_residual that stores references to the left hand side NiHu::function_space and the temporary NiHu::integral_transform. The additional functionality provided by the NiHu::weighted_residual class is that it can be evaluated into a matrix.

The next step is calling NiHu::operator<< with the matrix `W` as left hand side and the temporary NiHu::weighted_residual proxy as right hand side. NiHu::operator<< calls the member function template NiHu::weighted_residualNiHu::eval with the matrix as argument. This function is templated on the matrix type, so that weighted residuals can be evaluated into any indexable type (Eigen matrix of Matlab matrix, for example). The NiHu::weighted_residual::eval function simply forwards the result matrix reference and the stored test NiHu::function_space reference to the temporary NiHu::integral_transform's member function NiHu::integral_transform::test_on_into that _tests_ the integral transform `K[d]` _on_ the test function space `w` and writes the result _into_ the matrix `W`.

Function NiHu::integral_transform::test_on_into is able to determine during compilation if the stored trial function space `d` and the test function space parameter `w` are defined on the same mesh _type_. If no, then the underlying meshes are surely different and only regular integrals need to be evaluated. If yes, then singular integration may be needed. In the latter case, the mesh instances extracted from the function spaces are compared by address in order to make the decision.
\note This is why it is important to store function space views in references. If we make copies of meshes, then necessary singular integrals may be omitted during evaluation.

After this decision has been made, the appropriate instance of class template NiHu::assembly gets the scope. This class is responsible for assembling the output matrix `W` from element submatrices. Class NiHu::assembly is templated on the test NiHu::function_space type, the trial NiHu::function_space type and on the compile time boolean information whether singular integrals need to be computed or not. The static member function NiHu::assembly::eval_into is called with the NiHu::integral_operator `K`, the test function space `w`, the trial function space `d` and the result matrix `W` as arguments.

Let's take a break and summarise what happened up till now. We evaluated some operators, created proxy objects that were used to finally call NiHu::assembly::eval_into parametrised with the two function spaces `w` and `d`, the operator `K` and the result matrix `W`. The only runtime information processed was the comparison of two pointers (if the two mesh types are the same) in order to determine if we need to perform singular integrals or not. All remaining function calls and proxy instantiations are "optimised out" by the optimising compiler.

Evaluation on homogeneous subspaces {#tech_what_happens_dcall_each}
===================================

Function NiHu::assembly::eval_into splits the test and trial function spaces to homogeneous subspaces and matrix assembly is performed consecutively for each homogeneous subspace pair. For example: if the test function space consists of constant tria and constant quad fields and the trial function space consists of constant and isoparametric quads, then matrix assembly is performed separately for each row of the table below

Test space    | Trial space
--------------|--------------
constant tria | constant quad
constant tria | isoparametric quad
constant quad | constant quad
constant quad | isoparametric quad

This decomposition is performed during compile time, by utilising the library's code generation routine tmp::d_call_each. (d refers to a double nested loop over types)

The assembler's behaviour is slightly different if the integral operator `K` is a local operator (meaning that it acts only within an element, like the identity operator). In this case only those function space pairs are considered where the test and trial function spaces contain the same element type:

Test space    | Trial space
--------------|--------------
constant quad | constant quad
constant quad | isoparametric quad

Now the homogeneous assembler evaluates blocks of the result matrix `W` for each row of the tables, by traversing both homogeneous function spaces and evaluating the double integral

\f$
\displaystyle
W_{ij} = \left< w_i, \left(\mathcal{K}d_j\right)_{S_i}\right>_{F_j}
\f$

where \f$ S_i \f$ and \f$ F_j \f$ denote the subelements of the two (homogeneous sub)domains.

Evaluation of element integrals
===============================

The double integral in the above equation is evaluated by calling the function NiHu::integral_operator_base::eval_on_fields, which is specialised for the `Derived` class by the function NiHu::integral_operator::derived_eval_on_fields.
These specialisations take care of algebraic operations performed over integral operators, such as scalar multiplication.
In the special case of an identity operation, the evaluation is performed by the class NiHu::single_integral, whereas generally the same is achieved by using the class NiHu::double_integral.


Selecting the integral type
---------------------------

Depending on the formulation, single or double integrals are to be evaluated.
In case of double integrals the input is a test and a trial field and there are further cases that are treated separately.

In the general case (i.e. the kernel is singular and the test and trial spaces are interpreted over the same mesh) it is checked first if the integral to be computed is a singular or a regular one.
This check is performed by the evaluation of NiHu::element_match_eval function on the elements belonging to the test and trial fields.
The latter return an NiHu::element_match object, which stores the match type and the matching nodes if applicable.

1. If the integral is regular (i.e. the NiHu::element_match gives NiHu::NO_MATCH), the NiHu::complexity_estimator specified by the kernel class is evaluated in order to attain the quadrature degree needed for the integration.
When the quadrature degree is determined, the pre-evaluated shape function values are already available in the corresponding NiHu::store pattern.
This way shape functions values are evaluated only once for each regular quadrature, and this first and only evaluation happens when the corresponding _accelerator_ (NiHu::field_type_accelerator in this case) is first accessed.
(This pre-evaluation technique is referred to as _acceleration_ in the `NiHu` framework.)

2. If the integral is singular (i.e. the NiHu::element_match gives NiHu::FACE_MATCH, NiHu::EDGE_MATCH, or NiHu::CORNER_MATCH), the situation is bit more complicated.
Generally singular integrals can not be evaluated with good precision by means of standard quadratures, therefore their specialisations for different kernels, formalisms and element types is an essential part of the BEM.
To provide flexibility, two choices can be made when the kernel is defined.
	- When the kernel is **weakly singular** it means that the singularity can be cancelled by the appropriate choice of quadratures.
	The appropriate type of quadrature in this case is selectable by knowing only the _type of singularity_ (specified in the NiHu::singular_kernel_traits namespace) appearing in the kernel.
	Therefore this type of quadrature selection is referred to as **blind quadrature**, and is done automatically in `NiHu`.
	However, it is checked in compilation time if the specified kernel is integrable by means of a blind quadrature.
	If not, then the kernel is not weakly singular, the integral can not be performed, and the code will not compile.
	- To overcome the above limitations an alternative way is to specialise the NiHu::singular_integral_shortcut class for the given kernel, field types and match type, by which an integration method can be specified.
	If this specialisation of NiHu::singular_integral_shortcut exists, the blind quadrature branch above is omitted in compile time.
	This approach works for **weakly, strongly, and hypersingular** integrals.
	
A special case occurs when the kernel is not singular, or the integration is performed over two different meshes.
In this case the regular evaluation is automatic by skipping the NiHu::element_match_eval step.

Evaluation
----------

Both the regular and singular branches discussed above result in the evaluation of a block in the system matrix.
These blocks are calculated in most cases by iterating through a NiHu::quadrature, nevertheless, analytical integration in special cases can be implemented by means of the NiHu::singular_integral_shortcut.

When the double integral is performed using a quadrature rule, iterating through the quadratures over the test and trial fields is done by double cycle.
In the outer cycle the first kernel input (e.g. ::location) of the kernel is _bound_ to the input, that is calculated from the location of the current quadrature point on the physical element.
The result of this binding is a NiHu::kernel_base::kernel_bind type object, which can be evaluated given a single NiHu::kernel_input. (The NiHu::kernel_base::kernel_bind has an `operator()` function that takes only the second NiHu::kernel_input of the original kernel, while the first input is stored in the NiHu::kernel_base::kernel_bind object.)
Then, in the inner cycle this NiHu::kernel_base::kernel_bind object evaluates the kernel values, which are multiplied by the shape function values and quadrature weights and accumulated in the result matrix.
In `NiHu` the outer iteration is over the test field, while the inner cycle iterates through the quadrature over the trial field.

Matrix assembly
---------------

The accumulation process described above stores the values in a temporary matrix, which is allocated by the NiHu::double_integral::eval() function of the NiHu::double_integral class.
(The size of the temporary matrix are known from the number of test and trial field DOFs.)
The result matrix is returned to the NiHu::assembly class, which adds this matrix to a block of the overall system matrix, where the block indices are DOF indices.
The indexing is performed by the lightweight class NiHu::matrix_block.

The reason for evaluating into a smaller temporary matrix first is efficient memory management.
Since the overall system matrix can be huge (can occupy tens of gigabytes of memory), direct access to its elements corresponding to DOF indices that are relatively far can require paging, which could slow down the iteration of the accumulation process to a great extent.
By using the temporary matrix, the overall system matrix is accessed only once when an integral over a pair of fields is computed.
