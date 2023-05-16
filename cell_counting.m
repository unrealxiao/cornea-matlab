clear variables;
close all;
%address = '\\iu-opt-research\TankamLab\PPTs of RESULTS and PRESENTATIONS\POCM\Cell Counting\test_cell_2-1.tif';
address = 'G:\OCT_POCM_data\coopervision\Hsuanyeh_V1_2023-4-14\scan20\pic\reslice\AVG_crop.tif';
[I, cmap] = imread(address);
I = double(I);
%cali_factor = 1 / 0.001; %0.001 mm_meter per pixel, 1000 pixels per micro-meter
%I = double(rgb2gray(I));
I_r = I / max(I , [], 'all');
% [M,N]=size(I_r);
% I_r = I_r(N/2:N,:);
%Iblur = imbinarize(I_r, 0.4);
%Iblur = imgaussfilt(I_r, 1.5);
Iblur = I_r;%imlocalbrighten(I_r, 0.3);
%Iblur2 = imadjust(Iblur);
% % EM = graythresh(Iblur);
% [~, threshold] = edge(Iblur, 'sobel');
% 
% fudgeFactor = 0.5;
% BW = edge(Iblur,'sobel',threshold * fudgeFactor);
figure
imshow(Iblur, cmap)
%figure
%imshow(Iblur2)

%% Zero padding the image
%Nfft=1024;
Nfft = 1024;
max_side = max(size(Iblur));
fe= Nfft;

%Img = fft2(Iblur);
Img = fft2(Iblur, Nfft, Nfft);
figure
F = log(abs(fftshift(Img))).^4;
%F = abs(fftshift(Img));
%F = log(abs(Img));
%DF = fe/length(F); % frequency increment
% freqvec = -fe/2:DF:(fe/2);
freqvec = linspace(-fe/2,fe/2,Nfft);
imagesc(freqvec,freqvec,F)
Z = zeros(size(F));


F_dim = size(F);
center = round(F_dim/2);
%% construct each ring
% d = 2;%if I want 500 rings, what's the radii for each ring
% ring_x = center(1, 1):d:F_dim(1, 1); %the index of each rings in the F
% radius_ring = ring_x - center(1, 1);%the radius of each ring
% 
% radial_peak = zeros(2, length(radius_ring));%row 1 stores the adjusted frequency,
% %row 2 store the actual peak 
% radial_peak(1, :) = (radius_ring * 0.5) / center(1, 1);%adjust the radius to frequency 
% radial_peak = radial_peak(:, 1: end - 1);%remove the last column as there will be no data there
% %extract data from ring 
% 
% [Y,X]=ndgrid(1:F_dim(1, 2),1:F_dim(1, 1));
% delta_X = X-center(2);delta_Y = Y-center(1);%shift coordinate
% 
% for i = 1 : (length(ring_x) - 1)
%     inner_circle_radius = radius_ring(i);
%     outer_circle_radius = radius_ring(i + 1);
%     L = and(delta_X.^2 + delta_Y.^2 >= inner_circle_radius^2,  delta_X.^2 + delta_Y.^2 <= outer_circle_radius^2);%find index inside the ring(10, 20)
%     data_in_ring = F(L);
%     f = figure('visible','off');%hide the figure for later
%     H = histogram(data_in_ring);%histogram of data
%     [~, idx] = max(H.Values);%idx of the bin with highest density
%     Bin = H.BinEdges;%the bin in histogram
%     radial_peak(2, i) = Bin(idx);%this is the predominate frequency in the ring
% end
% 
% %% plot the radial graph 
% clear f H;
% figure;
% plot(radial_peak(1, :), radial_peak(2, :))
%imshow(L);%plot it, see if it  form a ring
%%
%radius = sqrt(center(1)^2 + center(2)^2);
radius = floor(fe / 2);
%radius = floor(Nfft) / 2;
dRho = 1;
%dRho = floor(Nfft / 2);
dTheta = 1 / radius;
%Thetas = (0:dTheta:2*pi);
Thetas = (0:dTheta:2*pi);
Rhos = (0:dRho:radius);
%Rhos = linspace(0,radius,dRho);

% polar mesh
[Theta, Rho] = meshgrid(Thetas, Rhos);

% transform...
[Xq,Yq] = pol2cart(Theta, Rho);

% translate to sit on the circle's center
Xq = Xq + center(2);
Yq = Yq + center(1);

% sample image at those points
radial_image = interp2(F, Xq, Yq);
figure;imagesc(radial_image)

%% radial average
% corp_row = min(center);
% corp_radial = radial_image(1:corp_row, :);%corp the longer side of the image
% imagesc(corp_radial)
%enhance_radial_image = imadjust(radial_image);
%average_radial = max(radial_image, [], 2);
average_radial = mean(radial_image, 2);
len_radial = size(average_radial);
x_radial = linspace(0, 0.5, len_radial(1, 1));
figure;
plot(x_radial(20:end), average_radial(20:end))
%plot(x_radial, average_radial)
%%
% crop_x = x_radial(100:end);
% crop_rad = fillmissing(average_radial(100:end), 'previous');
% f2 = fit(crop_x', crop_rad, 'exp1');
% plot(f2, crop_x, crop_rad)