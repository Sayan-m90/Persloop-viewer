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
