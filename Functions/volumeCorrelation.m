function corrVal = volumeCorrelation(vol1,vol2)

m1 = mean(vol1(:));
m2 = mean(vol2(:));

vol1mm = vol1(:)-m1;
vol2mm = vol2(:)-m2;

corrVal = sum(vol1mm.*vol2mm)/sqrt(sum(vol1mm.^2)*sum(vol2mm.^2));



end