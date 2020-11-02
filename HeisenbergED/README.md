**HeisenbergED** is an exact diagonalization implementation for the spin-1/2 XXX model, focusing on Sz- and momentum- sectors selection as a physics informed routine for the block-reduction of the eigenproblem. Due to prohibitive scaling of the many-body Hilbert space special attention has been devoted to memory management and fast build of matrices, following closely the standard reference by Sandvik ([arXiv:1101.3281](https://arxiv.org/abs/1101.3281)). Furthermore a statistical analysis based on *Phys. Rev. B 75, 155111 (2007)* has been carried out on the ED results. Correspondent preparatory lectures have been given by prof. Mario Collura at SISSA.

-----------------------

Material is partially coauthored with Francesca Paoletti <francesca.paoletti@sissa.it>, quoted here and there as *fp*.

More specifically: *She coded in f90. Here is reported only the data analysis, carried in python (+ some f90 commented code snippets, therein).
All the debugging has been carried out together, in a "cross-checked" fashion. Nevertheless the underlying logical structure of some subroutine is inherently distinct, due to authors' different taste.* 
