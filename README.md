# sqc-dft
ETO orbital based DFT with Lebedev quadrature (following the ADF paper)  - Simple quantum code (sqc) dft

For some project I needed to eval multicenter ETO integrals, I followed https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.1056 and eventually implemented many parts of the paper to have an actual DFT code.
ADF seems a very elegant and useful code, sadly it is very expensive to get the source code if you want to do experiment with your own ETO based theories. These implementations might help you.

Currently the only dependency is Eigen 3.3.5 (http://eigen.tuxfamily.org/index.php?title=Main_Page)

Later it will depend on libxc package to use realistic XC functional parametrizations.


Note: this project is not under active development. It is also lacking testing. Although many parts of the code, functions have been develepode alognside with Mathematica codes
 to ensure they are working.

The Lebedev quadrature is very elegantly takes care of the multicenter integrals, using only a few thousand points in the 3D space.

If you have interest in pushing this code forward or reuse it etc, then contact me, I am happy to help.

The idea is to be able to experiment with ETO style DFT code. ETO basis is intuitive in chemistry. Once an SFC Hamiltonian is obtained for a system,
 many physical quantities can be examined with non-interacting fashion, and if the SFC Hamiltonian is obtained with realistic coefficents the result is also expected to be realistic.

However keep in mind that for performant simulation this code is not suitable. (However it could be eventually)

TODO:

 - First of all, have a built in test case with solid results (for example compare with basic ADF results)
 - DIIS
 - better command line
 - write out outputs
 - libxc
 - Grid based openmp, not just atom, but need a table design..
 - Frozencore
 - Resolve issue with the same center integrals: if l,m != l,m we need to sum for grid points?
 - Mathematica based atomic test
 - http://webglmol.osdn.jp/glmol/embedding-examplesEN.html
