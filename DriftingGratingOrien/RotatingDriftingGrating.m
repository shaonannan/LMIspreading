%%% This matlab script is used for imshow population response for rotating
%%% drifting grating.
for i= 10:5:200
idx1 = 1:4:64;
figPha = reshape(mEbin_ra(i,idx1),4,4);
subplot(2,2,1);
imagesc(figPha);axis square;caxis([0 2.0]);
idx1 = 2:4:64;
figPha = reshape(mEbin_ra(i,idx1),4,4);
subplot(2,2,2);
imagesc(figPha);axis square;caxis([0 2.0]);
idx1 = 3:4:64;
figPha = figPha + reshape(mEbin_ra(i,idx1),4,4);
subplot(2,2,3);
imagesc(figPha);axis square;caxis([0 2.0]);
idx1 = 4:4:64;
figPha = figPha + reshape(mEbin_ra(i,idx1),4,4);
subplot(2,2,4);
imagesc(figPha);axis square;caxis([0 2.0]);
pause(0.1);
end