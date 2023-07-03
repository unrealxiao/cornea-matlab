clear variables
close all

%% model the surface of cornea and flatten it
%close all
[delta_matrix, smooth_surf, surface_cornea, Path_flat_save, Path_cross_flat, Path_cross_unflat, number_of_frames, original_scan, crop_scan, flip, Cl, cmap] = cornea_delta(25, 0.1, 0.6, "n");
%delta matrix return the distance between peak of the surface and each
%individual point
%(:, :, 1) of original_scan corresponds to first frame
%original_scan preserve the rows that were croped, while crop_scan deleted
%the row that were not interested. the reduced size in crop_scan help speed
%up the cross-corelation correction
% %% check the image
surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

figure;

surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% 
% disp('delta matrix complete')
% %% check the contrast at certain depth
% depth_frame = original_scan(:, 192, :);
% depth_frame = squeeze(depth_frame);
% imshow(depth_frame, cmap)
% [x, ~] = ginput(1);
% contrast = max(depth_frame, [], 'all') - min(depth_frame, [], 'all')
% if you want to flatten the data first, then use this code below

[flatten_original] = cornea_flatten(original_scan, delta_matrix, flip, Cl);%the flattened stack with row preserved
%[flatten_crop] = cornea_flatten(crop_scan, delta_matrix, flip, Cl);%the flattened stack with row croped
%use flatten_crop to do cross-correlation

%save flattened stack 

parfor i = 1 : number_of_frames
    sav_frame = flatten_original(:, :, i);
    
    if strcmp(Cl,'uint8')==1
        sav_frame=uint8(sav_frame); 
    elseif strcmp(Cl,'uint16')==1
        sav_frame=uint16(sav_frame); 
    end
    
    index = num2str(i);  % build filename with index
    imgname = strcat(Path_flat_save,'frameSh',index, '.DCM');

    dicomwrite(sav_frame,imgname);
    
end

disp('flattening complete')
%% let's see if substracting the top layer from the tissue layer will work

top_layer = squeeze(flatten_original(354, :, :));
tissue_layer = squeeze(flatten_original(366, :, :));
processed_layer = tissue_layer - top_layer;
%% we will then proceed to correct cross-correlation in the images and save the corrected frame

%we first correct the shift on the flattened images
[~, ~] = cross_correlation(flatten_crop, flatten_original, Path_cross_flat, flip, Cl);

disp('cross-correlation and saving complete')

%% then we can correct the non-flattened frame

[total_shift, shift_indice] = cross_correlation(crop_scan, original_scan, Path_cross_unflat, flip, Cl);

disp('cross-correlation and saving on unflattend complete')