function [total_shift, shift_indice] = cross_correlation(stack, full_stack, save_path, flip, Cl)

%stack should be croped size. This can speed up the cross-correlation
%process
size_stack = size(full_stack);
row_stack = size_stack(1, 1);
column_stack = size_stack(1, 2);
%correlation_correct_stack = zeros(size_stack);%creat another empty stack for correlation frame.

number_of_frames = size_stack(1, 3);

shift_indice = zeros(1, number_of_frames - 1); %store the shift for each frames

%we will determine the shift between all adjacent frames, and then store
%the shift in the shift_indice array. the shift_indice(1) will store the
%shift between flatten_stack(:, :, 1) and flatten_stack(:, :, 2)
stack_2 = stack; %create another flatten stack to avoid overhead in parloop
parfor j = 1 : (number_of_frames - 1)
    FF1 = stack(:, :, j);
    FF2 = stack_2(:, :, j + 1);
    l1=sum(FF1);
    %l2=sum(FF2);
    len=length(l1);
    
    %corel=xcorr(l1,l2);

    %figure; plot (corel);grid
    %[M, shift]=max(corel)

    Mcorel=xcorr2(FF1,FF2);
    [~, pos]=max(max(Mcorel))
    shift_indice(1, j)=pos-len;
end

disp('cross_correlation: correlation complete')
%compute the total shift for each frame
total_shift = zeros(1, number_of_frames - 1);

for n = 1 : (number_of_frames - 1)
    total_shift(1, n) = sum(shift_indice(1, 1 : n));%total_shift store the total shift for each frames
end
disp('cross_correlation: total_shift complete')
%correct correlation on the stack without row croped, and then save it

%save the first frame directly
first_frame = full_stack(:, :, 1);
if strcmp(Cl,'uint8')==1
    first_frame=uint8(first_frame); 
elseif strcmp(Cl,'uint16')==1
    first_frame=uint16(first_frame); 
end

if (flip == "Y")
   first_frame = flipud(first_frame);
end

index = num2str(1);  % build filename with index
imgname = strcat(save_path,'corl',index,'.DCM');



dicomwrite(first_frame,imgname);
%then save the rest of the frame

parfor m = 2 : number_of_frames
        
    corrected_frame = zeros(row_stack, column_stack);
    uncorrected_frame = full_stack(:, :, m);

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
        corrected_frame(:, 1 : (column_stack - pos_shift)) = uncorrected_frame(:, (pos_shift + 1) : column_stack);
        %corrected_frame(:, (column_stack - pos_shift + 1) : column_stack) = flatten_frame(:, 1 : pos_shift);
    else
        %corrected_frame(:, 1 : shift) = flatten_frame(:, (column_stack - shift + 1) : column_stack);
        corrected_frame(:, (shift + 1) : column_stack) = uncorrected_frame(:, 1 : (column_stack - shift));       
    end 
    
    %correlation_correct_stack(:, :, m) = corrected_frame;
    
    if strcmp(Cl,'uint8')==1
        corrected_frame=uint8(corrected_frame); 
    elseif strcmp(Cl,'uint16')==1
        corrected_frame=uint16(corrected_frame); 
    end
    
    if (flip == "Y")
       corrected_frame = flipud(corrected_frame);
    end
    
    index = num2str(m);  % build filename with index
    imgname = strcat(save_path,'corl',index,'.DCM');
    
    
    
    dicomwrite(corrected_frame,imgname);
    
end
end