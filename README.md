/*****************************************************************************/
"Persloop" software for computing persistent 1-cycle
/*****************************************************************************/

"Persloop" software is developed by the Jyamiti research group headed by 
Prof. Tamal K. Dey at the Department of Computer Science and Engineering 
of The Ohio State University.

The binaries are distributed for: Ubuntu Linux 64bit;

=================================
DESCRIPTION
=================================

Our software PersLoop: Given a barcode file, and a filtration file from the raw input, generates the persistent 1-cycle for the barcodes.
It implements Algorithm 3 of the following paper:
*****************************************************************************
Paper: T. K. Dey, T. Hao, S. Mandal
*Persistent 1-Cycles: Definition, Computation, and Its Application*
pg:123--136, Computational Topology in Image Context 2019, Springer International Publishing
*****************************************************************************
    
# Persloop-viewer
PersLoop is used to generate persistent 1-Cycles  
We address the problem of computing these representative cycles, termed as persistent 1-cycles, for H1-persistent homology with Z2 coefficients. 
We propose an alternative set of meaningful persistent 1-cycles that can be computed with an efficient polynomial time algorithm. 
We also inspect the stability issues of the optimal persistent 1-cycles and the persistent 1-cycles computed by our algorithm with the observation that the perturbations of both cannot be properly bounded. 
We design a software, presented in this webpage which applies our algorithm to various datasets. 
Experiments on 3D point clouds, mineral structures, and images show the effectiveness of our algorithm in practice.

This github containes all the related codes and softwares. 
Source code for Algorithm 3 shown in the paper is given in the folder perloop-src
persloop-bin contains precompiled Ubuntu librabries and persloop-example contains a file to run the script for its understanding
We suggest you generate barcodes using our previous software called simpers.(http://web.cse.ohio-state.edu/~dey.8/SimPers/Simpers.html)

The file Persloop-viewer-script contains the matlab script which gives correspondence between barcodes and the 1-cycle.
Meaning clicking on a barcode, you can see the corresponding barcode on the 3D off file or image.

If you want to work on dimensions~=3, compile and run using the folder persloop-src-all-dim. It will generate off files for the output cycles, but since these points are in higher dimension, you wont be able to visualise unless you use t-SNE or MDS or other dimension reduction software. You can use the coordinates, however, for higher dimension analysis.


The folder persloop-src contains source code,
The folder PersLoop-Viewer-Script contains an interactive Matlab program to view the cycles by clicking on the intervals(3D only)
