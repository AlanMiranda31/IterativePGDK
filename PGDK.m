clear all
close all
clc

addpath('Functions');
addpath('Datasets');

% Matlab .mat file containing variables "compositePrior" and "compositePriorReducedCounts"
priorsPETimageFileName = 'CompositePriorsExample_MLEM_mice_Fallypride.mat';

load(priorsPETimageFileName);

searchNeigh = 11; % closest neighbors cubic search window size. Must be odd
nknn = 100; % number of closest neighbors to consider for calculation of the kernel matrix
PGDmaxite = 1000; % maximum number of iterations for the PGD calculation
PGDtol = 0.00001; % convergence tolerance of the PGD optimization defined as the paramteres relative difference between iterations

kenelMatrixNS = PGDKernelMatrixCalc(compositePrior,compositePriorReducedCounts,nknn,searchNeigh,PGDmaxite,PGDtol);

% Denoised noisy (reduced counts) composite frames

filteredNoisyCompositeFrames = kenelMatrixNS*double(reshape(compositePriorReducedCounts,size(compositePriorReducedCounts,1)*size(compositePriorReducedCounts,2)*size(compositePriorReducedCounts,3),size(compositePriorReducedCounts,4)));
filteredNoisyCompositeFrames = reshape(filteredNoisyCompositeFrames,size(compositePrior,1),size(compositePrior,2),size(compositePrior,3),size(compositePrior,4));

