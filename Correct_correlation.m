clear variables

%% model the surface of cornea
%close all
[delta_matrix, smooth_surf, surface_cornea, Path_flat_save, Path_corl, number_of_frames, original_scan, crop_scan, flip, cmap, Cl] = cornea_delta(0.3, "n");
%delta matrix return the distance between peak of the surface and each
%individual point
%(:, :, 1) of original_scan corresponds to first frame
%original_scan preserve the rows that were croped, while crop_scan deleted
%the row that were not interested. the reduced size in crop_scan help speed
%up the cross-corelation correction

surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

 

figure;

surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')

disp('delta matrix complete')
%% if you want to flatten the data first, then use this code below

[flatten_original] = cornea_flatten(original_scan, delta_matrix, flip, Cl);%the flattened stack with row preserved
[flatten_crop] = cornea_flatten(crop_scan, delta_matrix, flip, Cl);%the flattened stack with row croped
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
%% we will then proceed to correct cross-correlation in the images 
size_stack = size(flatten_original); %
row_stack = size_stack(1, 1);
column_stack = size_stack(1, 2);
correlation_correct_stack = zeros(size_stack);%creat another empty stack for correlation frame.

[total_shift] = cross_correlation(flatten_crop);


disp('cross-correlation complete')

%% save the motion free frame and also the flattened non motion corrected frame
correlation_correct_stack(:, :, 1) = flatten_original(:, :, 1);

parfor m = 2 : number_of_frames
    corrected_frame = zeros(row_stack, column_stack);
    flatten_frame = flatten_original(:, :, m);
    
    shift = total_shift(m - 1);
    %assume positive shift mean the next frame move to the right
%     if shift > 0
%         corrected_frame(:, 1 : (column_stack - shift)) = flatten_frame(:, (shift + 1) : column_stack);
%         corrected_frame(:, (column_stack - shift + 1) : column_stack) = flatten_frame(:, 1 : shift);
%     else
%         positive_shift = abs(shift);
%         corrected_frame(:, 1 : positive_shift) = flatten_frame(:, (column_stack - positive_shift + 1) : column_stack);
%         corrected_frame(:, (positive_shift + 1) : column_stack) = flatten_frame(:, 1 : (column_stack - positive_shift));
%     end
    %assume positive shift move the next frame to the left
    if shift < 0
        pos_shift = abs(shift);
        corrected_frame(:, 1 : (column_stack - pos_shift)) = flatten_frame(:, (pos_shift + 1) : column_stack);
        %corrected_frame(:, (column_stack - pos_shift + 1) : column_stack) = flatten_frame(:, 1 : pos_shift);
    else
        %corrected_frame(:, 1 : shift) = flatten_frame(:, (column_stack - shift + 1) : column_stack);
        corrected_frame(:, (shift + 1) : column_stack) = flatten_frame(:, 1 : (column_stack - shift));       
    end
    
    %correlation_correct_stack(:, :, m) = corrected_frame;
    
    if strcmp(Cl,'uint8')==1
        corrected_frame=uint8(corrected_frame); 
    elseif strcmp(Cl,'uint16')==1
        corrected_frame=uint16(corrected_frame); 
    end
    
    index = num2str(m);  % build filename with index
    imgname = strcat(Path_corl,'corl',index,'.DCM');
    
    
    
    dicomwrite(corrected_frame,imgname);
    
end

disp('saving complete')

%% correct the non-flattened frame

