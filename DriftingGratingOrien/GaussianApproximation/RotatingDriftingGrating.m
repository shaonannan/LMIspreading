%%% This matlab script is used for imshow population response for rotating
%%% drifting grating.
for i= 40:5:2000
idx1 = 1:2:32;
figPha = reshape(mEbin_ra(i,idx1),4,4);
subplot(1,2,1);
imagesc(figPha);axis square;caxis([0 8.0]);
idx1 = 2:2:32;
figPha = reshape(mEbin_ra(i,idx1),4,4);
% idx1 = 3:2:32;
% figPha = figPha + reshape(mEbin_ra(i,idx1),4,4);
% idx1 = 4:2:32;
% figPha = figPha + reshape(mEbin_ra(i,idx1),4,4);
subplot(1,2,2);
imagesc(figPha);axis square;caxis([0 8.0]);
pause(0.1);
end