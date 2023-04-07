clear 
close all
%% test cornea_layer
%% crop images, decide the folder of images
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
dcm_file = dir([Path, '*.DCM']);
num_frames = length(dcm_file);
[I,cmap] = dicomread([Path, name]);

size_I = size(I);

imshow(I,cmap) 
lower_sd = 1.5;
upper_sd = 2;
row_frac = 4/6;
p = ginput(2);  % have user crop image by selecting 2 coordinate points 
p(p<1)=1;
% x1 = min(floor(p(1)), floor(p(2))); %xmin
y1 = min(floor(p(3)), floor(p(4))); %ymin
% x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
y2 = max(ceil(p(3)), ceil(p(4)));   %ymax
frameC2=I(y1:y2, :); 
figure
imshow(frameC2,cmap)
% num1 = 1;
% [miceye1, lower2, midd2, upper2] = cornea_layer(path, num1, 5, 3, 5);
 Bowman_3D = zeros(num_frames, size_I(1, 2));
%% lower_3D aquire
start_offset = 30;%start counting from the middle, instead of from the end to avoid 
%noise on the boundary of the images
bowman_ratio = 2;
parfor i = 1:num_frames
    [bowman_depth] = cornea_layer_notnormal(Path, i, y1, y2, start_offset, bowman_ratio);
    Bowman_3D(i, :) = bowman_depth;
end
%% only run this part once
[upper_row, upper_column, upper_v] = find(Bowman_3D);%extract all nonzero value and their indices in the matrix
%% flip the data(only run this once)
%upper_v = abs(size_I(1, 1) - upper_v);

%% fit the data by using griddata
[xq, yq] = meshgrid(1:num_frames, 1:size_I(1, 2));
upper_vq = griddata(upper_row, upper_column, upper_v, xq, yq);
%% eliminate all the data on the edge of the 2D data_matrix(so that the artifacts around the edge can be eliminated)

% for i = 1:1000
%     for k = 1:1000
%        if (i - 500)^2 + (k - 500)^2 > 300^2
%            lower_vq(i, k) = NaN;
%            midd_vq(i, k) = NaN;
%            upper_vq(i, k) = NaN;
%        end
%     end
% end

%% find the distance between epi and bowsman and then use smooth function
upper_vq = (y2 - y1 - upper_vq) + size_I(1, 1) - y2;
smooth_upper = smooth2a(upper_vq, 20, 20);

%% depth map
% epi_stroma_ko = smooth_upper - smooth_midd;
% epi_endo_ko = smooth_upper - smooth_lower;
% stroma_endo_ko = smooth_midd - smooth_lower;

%% plot the data

% mesh(epi_stroma, 'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','none')


surf(smooth_upper * 0.7,'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','none')

%zlim([200, 200])
%% plot the depth

 imagesc(smooth_upper * 0.7);title('epi');
 caxis([0 150]);

%%
real_epi_stroma = smooth_upper * 0.7;

%epi
thickness_epi = mean(real_epi_stroma, 'all', 'omitnan')%it seems like using mean(mean(A)) will give slightly different 
%result compared with mean(A, 'all')
std_epi = std(real_epi_stroma, 0, 'all', 'omitnan')
