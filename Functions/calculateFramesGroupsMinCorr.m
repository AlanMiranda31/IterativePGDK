function cpmM = calculateFramesGroupsMinCorr(Image4D,nframes,grs,ngrcalcite)

% Calculates groups of frames indices with minimum image correlation
% between the frames.
% 3D matrix "cpmM" contains in every row the indices of the corresponding
% group of frames, and is ordered along the 3rd dimension so that the first element ( i.e. cpmM(:,:,1) )has the
% least correlation among all frames groups, from all randomly created
% groups, and the last element ( i.e. cpmM(:,:,end) )
% the highest. To accomodate for groups with less than grs elements, the
% last row can have zeros

% ****** Input ****** 
% Image4D : dynamic PET image, in single precission 4D array format
% nframes : number of frames of the dynamic PET image
% grs : number of frames per group
% ngrcalcite : number of test groups 
%
% ****** Output ****** 
% cpmM :  3D matrix with frames groups indices

ngrp = ceil(nframes/grs);

framescorrM = zeros(nframes,nframes);

for i=1:nframes-1
    framei = Image4D(:,:,:,i);
    if sum(framei(:))>0
        for j=i+1:nframes
            framej = Image4D(:,:,:,j);
            corrVal = volumeCorrelation(framei,framej);
            framescorrM(i,j) = corrVal;
            
        end
    else
        for j=i+1:nframes
            framej = Image4D(:,:,:,j);
            if sum(framej(:))==0
                framescorrM(i,j) = 1;
            else
                framescorrM(i,j) = 0;
            end
        end
    end
end

framescorrM = framescorrM+framescorrM';
d = diag(framescorrM);
d(:) = 1;
framescorrM = framescorrM+diag(d);

% calculate mean correlation of frames groups

cpmM = zeros(ngrp,grs,ngrcalcite);
cmpC = zeros(ngrcalcite,1);

for i=1:ngrcalcite
    
    
    cpm = zeros(ngrp,grs);
    for k=1:ngrp
        ing = mod(k-1,ngrp)+1;
        if ing==1
            p = randperm(nframes);
        end
        
        pind = [(ing-1)*grs+1:ing*grs];
        pind(pind>nframes) = [];
        cp = p(pind);
        cpm(k,1:length(cp)) = cp;
    end
    
    zx = cpm(end,:);
    zx(zx==0) = [];
    if length(zx)==1
        cpm(end,2) = cpm(end-1,end);
        cpm(end-1,end) = 0;
    end
    
    cpmM(:,:,i) = cpm;
    
    grCorr = zeros(ngrp,1);
    
    for j=1:ngrp
        corrv = [];
        for k=1:grs-1
            f1 = cpm(j,k);
            if f1~=0
                for p=k+1:grs
                    f2 = cpm(j,p);
                    if f2~=0
                        corrv = [corrv framescorrM(f1,f2)];
                    end
                end
            end
        end
        
        
        grCorr(j) = mean(corrv);
        
    end
    
    cmpC(i) = mean(grCorr);

end

[cmpCsorted,cmpCsortedI] = sort(cmpC);

cpmM = cpmM(:,:,cmpCsortedI);

end