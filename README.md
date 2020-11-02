# CondMatLab @ [SISSA](https://cm.sissa.it)
All the scripts&codes for the 1st term courses I have taken at SISSA [2019/20 - PhD Program in Condensed Matter Theory].
Almost everything is coauthored, so no license is provided.

------------------------------------------------------------
Brief description of the contenent:

- *SolidStateLab* consists of three distinct projects, tackled within the infamous *Solid-State Problems* course given by prof. Stefano de Gironcoli at SISSA. Main focus is electronic structure theory, ranging from more-or-less-realistic tight binding calculations on 2D materials to berry-phase-related properties of minimalistic topological models. The final problem set on lattice dynamics shifts towards a more quantitative point of view, focusing on accurate reproduction of experimental literature.

- *HFHubbardLab* is a theoretically inclined training on static mean-field (*aka* Hartree Fock) treatment of Hubbard models. Nevertheless some interesting numerical aspects arise as stability of the self-consistency loop is concerned (mainly around first-order critical lines, where the loop may converge on metastable solutions). The correspondent preparatory course has been given by prof. Massimo Capone at SISSA.

- *HeisenbergED* is an exact diagonalization implementation for the spin-1/2 XXX model, focusing on Sz- and momentum- sectors selection as a physics informed routine for the block-reduction of the eigenproblem. Due to prohibitive scaling of the many-body Hilbert space special attention has been devoted to memory management and fast build of matrices, following closely the standard reference by Sandvik ([arXiv:1101.3281](https://arxiv.org/abs/1101.3281)). Furthermore a statistical analysis based on *Phys. Rev. B 75, 155111 (2007)* has been carried out on the ED results. Correspondent preparatory lectures have been given by prof. Mario Collura at SISSA.

Side note:

- [*dMFTlab*](https://github.com/Bellomia/dMFTlab) could have been part of this repository, for it has started as a "snooping around" stage within the Advanced Many-Body course given by Massimo Capone. But I decided to host it elsewhere as an active and growing playground to give research ideas a quick and inexpensive test.

