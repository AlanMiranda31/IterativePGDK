clear all
close all
clc

addpath('Functions');
addpath('Datasets');

% Matlab .mat file containing varaible "Image4D" 
dynamicPETimageFileName = 'Image4Dexample_MLEM_mice_Fallypride.mat';

load(dynamicPETimageFileName);

searchNeigh = 11; % closest neighbors cubic search window size. Must be odd
nknn = 100; % number of closest neighbors to consider for calculation of the kernel matrix
framesGroupSize = 4; % frames group size for calculation of the iterative PGD
numberofItePGDKiterations = 10; % total number of iterations for the iterative PGD algorithm
PGDmaxite = 1000; % maximum number of iterations for the PGD calculation
PGDtol = 0.00001; % convergence tolerance of the PGD optimization defined as the paramteres relative difference between iterations

[itePGDkenelMatrix,Image4Dfiltered] = iterativePGDKernelMatrixCalcProj(Image4D,framesGroupSize,numberofItePGDKiterations,nknn,searchNeigh,PGDmaxite,PGDtol);

