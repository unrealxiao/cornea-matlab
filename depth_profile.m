clear variables;
close all;
address_a = 'D:\xiaoliu_onedrive\OneDrive - Indiana University\lab\POCM-polarization measure\power supply depth profile\';
name_1 = 'source_burn_r45';
name_2 = 'source_burn_r0';
name_3 = 'source_burn_r225';
fileID = fopen([address_a, name_1], 'r');
fileID2 = fopen([address_a, name_2], 'r');
fileID3 = fopen([address_a, name_3], 'r');
formatSpec = '%f';
A = fscanf(fileID, formatSpec);
B = fscanf(fileID2, formatSpec);
C = fscanf(fileID3, formatSpec);
A = A(1:11:end);
B = B(1:11:end);
C = C(1:11:end);
x_scale = (1:length(A)/2) * 0.00096966;% scaling factor in air
% [image1, cmap1] = imread([address_a, strcat(name_1, num2str(name_2), name_3)]);
% [image2, cmap2] = dicomread([address_a, strcat(name_1, num2str(name_2 - 1), name_3)]);
% [image3, cmap3] = dicomread([address_a, strcat(name_1, num2str(name_2), name_3)]);
% [image4, cmap4] = dicomread([address_a, strcat(name_1, num2str(name_2 + 1), name_3)]);
% [image5, cmap5] = dicomread([address_a, strcat(name_1, num2str(name_2 + 2), name_3)]);

%% plot image and its depth map
savefolder = '\\iu-opt-research\TankamLab\PPTs of RESULTS and PRESENTATIONS\POCM\Calibration\cross-correlation-images';

close all
figure;
plot(x_scale, A(2049:end), "blue")
title('0 degree')
xlabel('Depth(mm)', 'fontweight','bold')
ylabel('Intensity(a.u)', 'fontweight','bold')
ylim([0 10]);
ax = gca;
ax.FontWeight = 'bold';
saveas(gcf, [savefolder, '\0degree.png'])

figure;
plot(x_scale, B(2049:end),  "green")
title('45 degree')
xlabel('Depth(mm)', 'fontweight','bold')
ylabel('Intensity(a.u)', 'fontweight','bold')
ax = gca;
ax.FontWeight = 'bold';
ylim([0 2.5]);
xlim([0.5 0.58]);
saveas(gcf, [savefolder, '\45degree.png'])

figure;
plot(x_scale, C(2049:end),  "red")
title('22.5 degree')
xlabel('Depth(mm)', 'fontweight','bold')
ylabel('Intensity(a.u)', 'fontweight','bold')
ax = gca;
ax.FontWeight = 'bold';
ylim([0 10]);
saveas(gcf, [savefolder, '\225degree.png'])
%contrast = max(Avg_img, [], 'all') - min(Avg_img, [], 'all');

%% zoom
close all
figure;
half_axis = length(x_scale);
plot(x_scale, A(2049:end), "blue")
hold on
plot(x_scale, B(2049:end),  "green")
hold on
plot(x_scale, C(2049:end),  "red")
xlabel('Depth(mm)', 'fontweight','bold')
ylabel('Intensity(a.u)', 'fontweight','bold')
ax = gca;
ax.FontWeight = 'bold';
ylim([0 7]);

xlim([0.3 0.6]);
%xlim([0.16 0.36])


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


