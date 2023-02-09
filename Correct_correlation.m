clear variables

%% model the surface of cornea
%close all
[delta_matrix, smooth_surf, surface_cornea, Path_flat_save, Path_corl, number_of_frames, original_scan, flip, cmap, Cl] = Cornea_flatten(0.3, "Y");
%delta matrix return the distance between peak of the surface and each
%individual point
%(:, :, 1) of original_scan corresponds to first frame
surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

 

figure;

surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')

cat('delta matrix complete')
%% if you want to flatten the data first, then use this code below
size_stack = size(original_scan);
row_stack = size_stack(1, 1);
column_stack = size_stack(1, 2);
flatten_stack = zeros(size_stack);%create another empty stack for flattened frame

%flatten the image, and put flattened frame into the 3 dimensional array

parfor k = 1:number_of_frames
    
    original_frame = original_scan(:, :, k);
    frameShtwo=zeros(row_stack, column_stack); 
    
    %moving pixel on original_frame
        
    for col = 1:column_stack
        if delta_matrix(k, col) > 0
            frameShtwo(1:(row_stack - delta_matrix(k, col)), col) = original_frame((delta_matrix(k, col) + 1):row_stack, col);
            frameShtwo((row_stack - delta_matrix(k, col) + 1):row_stack, col) = original_frame(1:delta_matrix(k, col), col);
        else
            positive_delta = abs(delta_matrix(k, col));
            frameShtwo(1:positive_delta, col) = original_frame((row_stack - positive_delta + 1):row_stack, col);
            frameShtwo((positive_delta + 1):row_stack, col) = original_frame(1:(row_stack - positive_delta), col);
        end 
    end
    
  
    
    
    
    if strcmp(Cl,'uint8')==1
        frameShtwo=uint8(frameShtwo); 
    elseif strcmp(Cl,'uint16')==1
        frameShtwo=uint16(frameShtwo); 
    end

    

    %frameShtwo = imgaussfilt(frameShtwo, 0.7);%use gaussian filter to smooth the small discontinuities caused by flattening.
    %frameShtwo = medfilt2(frameShtwo, [2, 2]);%using moving average
    if (flip == "Y")
       frameShtwo = flipud(frameShtwo);
    end
    
    flatten_stack(:, :, k) = frameShtwo;%put flattened frame into flattened stack

  
end

cat('flattening complete')
%% we will then proceed to correct the horizontal motion in the images 

correlation_correct_stack = zeros(size_stack);%creat another empty stack for correlation frame.

shift_indice = zeros(number_of_frames - 1); %store the shift for each frames

%we will determine the shift between all adjacent frames, and then store
%the shift in the shift_indice array. the shift_indice(1) will store the
%shift between flatten_stack(:, :, 1) and flatten_stack(:, :, 2)
flatten_stack_2 = flatten_stack; %create another flatten stack to avoid overhead in parloop
parfor j = 1 : (number_of_frames - 1)
    FF1 = flatten_stack(:, :, j);
    FF2 = flatten_stack_2(:, :, j + 1);
    l1=sum(FF1);
    l2=sum(FF2);
    len=length(l1);
    
    corel=xcorr(l1,l2);

    %figure; plot (corel);grid
    %[M, shift]=max(corel)

    Mcorel=xcorr2(FF1,FF2);
    [m, pos]=max(max(Mcorel))
    shift_indice(j)=pos-len;
end

%compute the total shift for each frame
total_shift = zeros(number_of_frames - 1);

for n = 1 : length(shift_indice)
    total_shift(n) = sum(shift_indice(1 : n));%total_shift store the total shift for each frames
end
%correct the shift for each frame

cat('cross-correlation complete')

%% save the motion free frame
correlation_correct_stack(:, :, 1) = flatten_stack(:, :, 1);

parfor m = 2 : number_of_frames
    corrected_frame = zeros(row_stack, column_stack);
    flatten_frame = flatten_stack(:, :, m);
    
    shift = total_shift(m - 1);
    
    if shift > 0
        corrected_frame(:, 1 : (column_stack - shift)) = flatten_frame(:, (shift + 1) : column_stack);
        corrected_frame(:, (column_stack - shift + 1) : column_stack) = flatten_frame(:, 1 : shift);
    else
        positive_shift = abs(shift);
        corrected_frame(:, 1 : positive_shift) = flatten_frame(:, (column_stack - positive_shift + 1) : column_stack);
        corrected_frame(:, (positive_shift + 1) : column_stack) = flatten_frame(:, 1 : (column_stack - positive_shift));
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

cat('saving complete')

%% correct the non-flattened frame

