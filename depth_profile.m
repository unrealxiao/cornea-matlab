address_a = 'G:\OCT_POCM_data\Glass-for-double-peak-gaps-measure\glass-R45\pic\';
name_1 = 'frame';
name_2 = 250;
name_3 = '.DCM';
[image1, cmap1] = dicomread([address_a, strcat(name_1, num2str(name_2 - 2), name_3)]);
[image2, cmap2] = dicomread([address_a, strcat(name_1, num2str(name_2 - 1), name_3)]);
[image3, cmap3] = dicomread([address_a, strcat(name_1, num2str(name_2), name_3)]);
[image4, cmap4] = dicomread([address_a, strcat(name_1, num2str(name_2 + 1), name_3)]);
[image5, cmap5] = dicomread([address_a, strcat(name_1, num2str(name_2 + 2), name_3)]);

Avg_img = (image1 + image2 + image3 + image4 + image5) / 5;
[row, col] = size(Avg_img);
%figure;
plot(Avg_img(:, round(col/2)))
ylim([0 250])
hold on
%% plot image and its depth map
figure;
imshow(Avg_img, cmap1);


%% find the difference, sum, and ratio of I0 and I45
address_45 = 'G:\OCM_FM\E\3Dimage\coopervision\Exp_Sep_16_2022\IRCard\Scan1_S45_R45\pic1\';
address_0 = 'G:\OCM_FM\E\3Dimage\coopervision\Exp_Sep_16_2022\IRCard\Scan1_S45_R0\pic1\';

folder_sum = 'G:\OCM_FM\E\3Dimage\coopervision\Exp_Sep_16_2022\IRCard\sum_scan1'; %store I45 + I0
folder_diff = 'G:\OCM_FM\E\3Dimage\coopervision\Exp_Sep_16_2022\IRCard\diff_scan1';%store I45 - I0
folder_ratio = 'G:\OCM_FM\E\3Dimage\coopervision\Exp_Sep_16_2022\IRCard\ratio_scan1';%store I45 - I0 / I45 + I0
if ~exist(folder_sum, 'dir')
    mkdir(folder_sum)
end

if ~exist(folder_diff, 'dir')
    mkdir(folder_diff)
end

if ~exist(folder_ratio, 'dir')
    mkdir(folder_ratio)
end

name_1 = 'frame';
name_2 = 500;
name_3 = '.DCM';
[test_image, cmap3] = dicomread([address_45, strcat(name_1, num2str(name_2), name_3)]);
Cl=class(test_image);

parfor i = 1:1000
   index = num2str(i);  % build filename with index
   filename = strcat(name_1,index); 
   [frame_45,~] = dicomread([address_45, filename]); % open frame
   [frame_0,~] = dicomread([address_0, filename]);
   frame_sum = frame_45 + frame_0;
   frame_diff = frame_45 - frame_0;
   frame_ratio = rdivide(frame_diff,frame_sum);
    if strcmp(Cl,'uint8')==1
        frame_sum = uint8(frame_sum); 
        frame_diff = uint8(frame_diff); 
        frame_ratio = uint8(frame_ratio); 
    elseif strcmp(Cl,'uint16')==1
        frame_sum = uint16(frame_sum); 
        frame_diff = uint16(frame_diff); 
        frame_ratio = uint16(frame_ratio); 
    end  
    
    sum_save = [folder_sum, '\', name_1, index, '.DCM'];
    diff_save = [folder_diff, '\', name_1, index, '.DCM'];
    ratio_save = [folder_ratio, '\', name_1, index, '.DCM'];
    
    dicomwrite(frame_sum, sum_save);
    dicomwrite(frame_diff, diff_save);
    dicomwrite(frame_ratio, ratio_save);
end


