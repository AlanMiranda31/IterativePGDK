
Author: Alan Miranda 
alan.mirandamenchaca@uantwerpen.be
University of Antwerp     
2025

Code to perform calculation of the kernel matrix for (3D) PET dynamic reconstruction
using the projected gradient descent (PGD) and the iterative PGD methods.

%%%%%%%%%%%%%%%%%%%%%%

PGD kernel matrix

The inputs to the algorithm are:
- High quality PET composite images, defined as 4D images where the 3D composite images are concatenated along the 4th dimension
- Reduced counts PET composite images, corresponding to the respective high quality composite image
- Nearest neighbors search window, defining a square search neighborhood
- Number of nearest neighbors
- Maximum number of iterations for the PGD optimization
- Convergence tolerance of the PGD optimization, defined as the relative difference in parameters change between iterations

The output is:
- The kernel matrix, provided in a Matlab sparse matrix of size NVOX x NVOX, where NVOX is the number of voxels in the 3D PET image

%%%%%%%%%%%%%%%%%%%%%%

Iterative PGD kernel matrix

The supplied code is for the method using the (simple) kernel matrix multiplication approach to denoise the PET image frames at every iteration. If
you want to follow the approach using the KEM reconstruction denoising, you can insert your code for KEM reconstruction at the end of the main for loop.

The inputs to the algorithm are:
- The (noisy) independent frame 4D dynamic PET reconstruction with the (final) target framing, where frames are concatenated in the 4th dimension
- Number of frames per group to calculate the PGD kernel matrix
- Number of iterations for the iterative PGD kernel matrix calculation
- Nearest neighbors search window, defining an square search neighborhood
- Number of nearest neighbors
- Maximum number of iterations for the PGD optimization
- Convergence tolerance of the PGD optimization, defined as the relative difference in parameters change between iterations

The output is:
- The kernel matrix, provided in a Matlab sparse matrix of size NVOX x NVOX, where NVOX is the number of voxels in the 3D PET image
- The filtered dynamic PET image calculated in the last iteration as the product of the kernel matrix and the dynamic PET image frames

%%%%%%%%%%%%%%%%%%%%%%

Datasets examples are supplied in Matlab '.mat' file format, corresponding to a dynamic PET [18F]Fallypride scan of 2 mice simultaneoulsy.

File "CompositePriorsExample_MLEM_mice_Fallypride.mat" contains the composite priors needed for the PGD kernel matrix calculation, while
"Image4Dexample_MLEM_mice_Fallypride.mat" contain the MLEM dynamic PET image of the same scan, needed for the iterative PGD kernel matrix calculation.

Precompiled mex files, for windows and linux, are supplied for the CPU parallel calculation of the nearest neighbors and for the PGD optimization.






