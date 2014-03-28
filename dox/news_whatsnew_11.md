What's new in release 1.1 {#news_whatsnew_11}
=========================

\page news_whatsnew_11

[guiggi]: \ref tut_guiggiani_hypersingular "Hypersingular integrals and Guiggiani's method"
[vector]: \ref tut_funspace_vector "Vector function spaces"
[customkernel]: \ref tut_custom_kernel "Custom kernel definitions"

**Function space representations**
- Vector fields and field generation options (see [vector])
- Automatic store or on the fly decisions
- Traits namespaces instead of traits classes
- Static members in libraries (lib)

**Kernels**
- Support for easily defining single-brick kernels with a functor object (see [customkernel])
- Asymptotic behaviour defined
- Library: elastostatics
- Library: interval estimators for log type singularities

**Integration and Matrix assembly**
- Guiggiani's method (see [guiggi])
- Library: Specialised for scalar potential problems and elastostatics
- Library: Acoustic Burton-Miller approach with non-planar elements
- Compile time check of possible match types
- Block product: vector problems 

