function [itePGDkenelMatrix,Image4Dfiltered] = iterativePGDKernelMatrixCalcProj(Image4D,grs,iteknite,nknn,searchNeigh,PGDmaxite,PGDtol)

% Calculates the (sparse) kernel matrix for PET dynamic reconstructions
% with the iterative PGD algorithm (itePGDK)
%
% The liner voxel index is calculated following Matlab indexing, i.e.
% column-wise. For example the linear voxel index of image element at row = 13,
% column = 11, depth = 21, in an image with xsize = 128, ysize = 128, and zsize = 159 is:
%
% linear_index = row+(ysize)*(column-1)+(xsize*ysize)*(depth-1) = 328973
%
% Sparse kernel matrix row 328973 therefore has the kernel matrix weights of this
% voxel
%
% ****** Input ****** 
% Image4D : Dynamic PET image, in single precission 4D array format, calculated with e.g. MLEM reconstruction
% grs : frames group size for the itePGDK calculation
% iteknite : total number of iterations to run the itePGDK calculation
% nknn : number of closest neighbors for the kernel matrix calculation
% searchNeigh : closest neighbors cubic window size. Must be odd
% PGDmaxite : maximum number of iterations for the PGD calculation
% PGDtol : convergence tolerance of the PGD optimization defined as the paramteres relative difference between iterations
%
% ****** Outut ****** 
% itePGDkenelMatrix : kernel matrix in sparse Matlab matrix format
% Image4Dfiltered : the filtered Image4D by multiplication with itePGDkenelMatrix in the last iteration

nframes = size(Image4D,4);

xsize = size(Image4D,2);
ysize = size(Image4D,1);
zsize = size(Image4D,3);

nvox = xsize*ysize*zsize;

meanImage = zeros(ysize,xsize,zsize,'single');
meanImage(:,:,:,2:end) = [];
for i=1:size(Image4D,4)
    meanImage = meanImage+Image4D(:,:,:,i);
end

disp('Selecting non-background voxels')
% Image segmentation based on Otsu threshold to obtain active voxels

segmentationthresh = 0.1;
otth = calcOtsuThresholdSingle(meanImage);
mask = zeros(size(meanImage),'single');
mask(meanImage(:)>otth*segmentationthresh) = 1;
CC = bwconncomp(mask);
maskTemp = zeros(size(mask),'single');
for i=1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{i})>1000
        maskTemp(CC.PixelIdxList{i}) = 1;
    end
end
maskFilled = imfill(maskTemp,'holes');


inNoAct = find(maskFilled(:)==0);

Image4Dfiltered = Image4D;

featdyn = reshape(Image4D,nvox,nframes)';
featdynFilt = reshape(Image4Dfiltered,nvox,nframes)';

flagInAct = ones(nvox,1,'int32');
flagInAct(maskFilled(:)==0) = 0;

% calculation of frames groups with minimum image correlation between frames
% within the group
ntestgroups = 10000;
cpmM = calculateFramesGroupsMinCorr(Image4D,nframes,grs,ntestgroups);

nvoxact = sum(flagInAct);
indvoxact = find(flagInAct==1);

%% run kernel matrix calculation
disp('Starting calculation of itePGDK kernel matrix')
for i=1:iteknite
    disp(['Running iteration ',num2str(i),' of ',num2str(iteknite) ,' of the iterative PGD kernel matrix calculation']);
    cpm = cpmM(:,:,i);
  
    % calculation of the nearest neighbors in groups followed by
    % calculation of the most repited indices
    Idxknn = findKNNgroups(featdynFilt,int32(xsize),int32(ysize),int32(zsize),int32(searchNeigh),int32(nknn),int32(maskFilled),int32(cpm));
    
    kernelMatrixInit = ones(size(Idxknn),'single');
    kernelMatrixInit(:) = 1/nknn;

    kenelMatrixAdjustedSum = zeros(size(kernelMatrixInit),'single');
    
    % calculation of the kernel matrix for all groups
    for jj=1:size(cpm,1)
        
        cp = cpm(jj,:);
        cp(cp==0) = [];
        
        featdynSmall = featdyn(cp,:);
        featdynFiltSmall = featdynFilt(cp,:);
        
        [kenelMatrixAdjusted,objFvalAccC] = kernelPGDcalculation(featdynFiltSmall,featdynSmall,Idxknn,...
            kernelMatrixInit,int32(PGDmaxite),single(PGDtol),flagInAct);
        
        for j=1:nvox
            if flagInAct(j)==1
                n = sum(kenelMatrixAdjusted(j,:));
                if n~=0
                    kenelMatrixAdjusted(j,:) = kenelMatrixAdjusted(j,:)/n;
                end
            end
        end
        
        kenelMatrixAdjustedSum = kenelMatrixAdjustedSum+kenelMatrixAdjusted;

    end
    
    % average of all groups kernel matrix
    kenelMatrixAdjusted = kenelMatrixAdjustedSum/size(cpm,1);
    
    for j=1:nvox
        if flagInAct(j)==1
            n = sum(kenelMatrixAdjusted(j,:));
            if n~=0
                kenelMatrixAdjusted(j,:) = kenelMatrixAdjusted(j,:)/n;
            end
        end
    end
    
    % storing of the most repeated nearest neighbors and correspondig
    % matrix weights for every iteration
    
    if i==1
        kernelKNNind = cell(nvoxact,1);
        kernelKNNindN = cell(nvoxact,1);
        kernelKNNw = cell(nvoxact,1);
        
        for j=1:nvoxact
            kernelKNNind{j} = Idxknn(indvoxact(j),:);
            kernelKNNindN{j}  = ones(1,nknn);
            kernelKNNw{j}  = kenelMatrixAdjusted(indvoxact(j),:);
        end
    else
        for j=1:nvoxact
            
            vind = kernelKNNind{j};
            nind = kernelKNNindN{j};
            wind = kernelKNNw{j};
            
            
            iind = Idxknn(indvoxact(j),:);
            iw = kenelMatrixAdjusted(indvoxact(j),:);
            
            [C,ia,ib] = intersect(vind,iind);
            
            if ~isempty(C)
                nind(ia) = nind(ia)+1;
                wind(ia) = wind(ia)+kenelMatrixAdjusted(indvoxact(j),ib);
            end
            
            [C,ib] = setdiff(iind,vind);
            
            if ~isempty(C)
                vind = [vind C];
                nind = [nind ones(1,length(C))];
                wind = [wind kenelMatrixAdjusted(indvoxact(j),ib)];
            end
            
            kernelKNNind{j} = vind;
            kernelKNNindN{j} = nind;
            kernelKNNw{j} = wind;
            
            
        end
        
    end
    
    % calculation of the kernel matrix weights based on most repeated
    % nearest neighbors indices
    
    Idxknn = zeros(nvox,nknn,'int32');
    Idxknn(:,1) = [1:nvox]';
    
    for j=1:nvoxact
        
        vind = kernelKNNind{j};
        nind = kernelKNNindN{j};
        wind = kernelKNNw{j};
        
        
        wind = wind./single(nind);
        wind(isnan(wind)) = 1/nknn;
        
        [nindSorted,nindSortedI] = sort(nind,'descend');
        
        Idxknn(indvoxact(j),:) = vind(nindSortedI(1:nknn));
        kenelMatrixAdjusted(indvoxact(j),:) = wind(nindSortedI(1:nknn));
        
    end
    
    
    for j=1:nvox
        if flagInAct(j)==1
            n = sum(kenelMatrixAdjusted(j,:));
            if n~=0
                kenelMatrixAdjusted(j,:) = kenelMatrixAdjusted(j,:)/n;
            end
        end
    end
    
    kernelMatrix = kenelMatrixAdjusted;
    kernelMatrix(inNoAct,1) = 1.0;
    kernelMatrix(inNoAct,2:end) = 0.0;
    
    zz = find(kernelMatrix(:)>0);
    IdxknnL = Idxknn(zz);
    kernelMatrixL = kernelMatrix(zz);
    auxrow = repmat([1:nvox]',1,nknn);
    auxrow = auxrow(zz);
    itePGDkenelMatrix = sparse(auxrow(:),double(IdxknnL(:)),double(kernelMatrixL(:)),nvox,nvox);
    
    featdynFilt = itePGDkenelMatrix*double(featdyn');
    featdynFilt = single(featdynFilt');
    
    % check which voxels have too many kernel weights close to zero and
    % mark them in the "maskFilled" image
    
    % threshold corresponding to the percentage of the maximum to consider
    % a small weight kernel value
    bfth = 0.05;
    
    for j=1:nvoxact
        
        kenelMatrixAdjustedTemp  = kenelMatrixAdjusted(indvoxact(j),:);
        
        kenelMatrixAdjustedTemp1 = ones(size(kenelMatrixAdjustedTemp));
        
        [vmax,imax] = max(kenelMatrixAdjustedTemp);
        
        kenelMatrixAdjustedTemp1(kenelMatrixAdjustedTemp<vmax*bfth) = 0;
        
        if sum(kenelMatrixAdjustedTemp1)<nknn/2%indvoxact(j)==selVox1
            maskFilled(indvoxact(j)) = maskFilled(indvoxact(j))+1;
        end
        
    end
    
    %%%%%% To use the KEM denosing approach add here your code to calculate
    %%%%%% the KEM reconstruction. The result of this calculation will
    %%%%%% replace the variable featdynFilt
    
end
disp('Calculation of the itePGDK kernel matrix finished')
Image4Dfiltered = single(reshape(featdynFilt',ysize,xsize,zsize,nframes));


end