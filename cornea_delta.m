function [delta_matrix, smooth_surf, smooth_poly, Path_save_flat, Path_save_cross_flat, Path_save_cross_unflat, numframes, original_stack, crop_stack, flip, Cl, cmap] = cornea_delta(start_offset, top_ratio, low_ratio, flip)
%if you want to flip the images(so that endothelium will be on top), enter
%flip as "Y"
%top_ratio represents the top percents of the column array, lower_ratio
%represents the low percents of the column array. use to find the peak
%start offset use to determine the gap between the peak point and the
%array we are going to average. This is used to make sure that the array is
%not on the hill that leads to the peak. without this value the peak won't
%be significanly higher than the average, which leads to the failure of
%identifying the peak.

%smooth_surf are the surface of the cornea, row 1 represents the segment of
%frame1, row2 to frame2, and so on
Path=[pwd,'/'];
[~,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
%Ext=name(end-3:end);
mkdir(Path,'Flat'); % Create a subfolder in the directeory
%mkdir(Path, 'cross_flat');
%mkdir(Path, 'cross_unflat')
Path_save_flat=[Path,'Flat\'];
Path_save_cross_flat = [Path, 'cross_flat\'];
Path_save_cross_unflat = [Path, 'cross_unflat\'];
%count the number of files in the directories
dcm_file = dir([Path, '*.DCM']);

numframes = length(dcm_file); % number of input frames

num=1; %initial frame number 
%deg=3; %degree of the polynomial fit
%deg2=2;
%threshold=thres_hold; % for the peak detection 
name=[Path,append('frame', num2str(round(numframes/2)), '.DCM')];
[I,cmap] = dicomread(name);  % find out size of images by importing one
if (flip == "Y")
   I = flipud(I); 
end
%I=flipud(I);


Cl=class(I);
tic

%% crop

frame = I;

size_of_frame = size(frame);

imshow(frame,cmap) 
p = ginput(2);  % have user crop image by selecting 2 coordinate points 
p(p<1)=1;
%x1 = min(floor(p(1)), floor(p(2))); %xmin
y1 = min(floor(p(3)), floor(p(4))); %ymin
%x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
y2 = max(ceil(p(3)), ceil(p(4)));   %ymax

%frameC2=frame(y1:y2, x1:x2);
frameC2=frame(y1:y2, :);
imshow(frameC2,cmap) %frameC = cropped frame 

size_of_crop = size(frameC2);

S=size(frameC2);
%H=S(1,1); %height of figure (rows) 
L=S(1,2); %lenght of figure (columns) 

%%
name_5= 'frame';

% Perform the parallel computing

polysf=zeros(numframes,L); %yy4 surface reconstruction  

%create 3 dimensional array to store the whole stacks

original_stack = zeros(size_of_frame(1), size_of_crop(2), numframes);%store the frames that preserve row

crop_stack = zeros(size_of_crop(1), size_of_crop(2), numframes);%use to compute the cross-corelation later

%flattened_stack = zero(size_of_frame(1, 1), size_of_frame(1, 2), numframes);

parfor i = 1:numframes

    count=i+num-1;
    %count = 300;

    index = num2str(count);  % build filename with index
    filename = strcat(name_5,index);
    
    [frame,~] = dicomread([Path, filename]);% open frame
    %frame = imbinarize(frame, 'adaptive', 'ForegroundPolarity', 'bright', 'Sensitivity', 0.2);% binarization
    %frame = medfilt2(frame);
    %frame = rot90(frame,0); % rotate frames, this can be used if the frames are not in the proper orientation
    if (flip == "Y")
       frame=flipud(frame); 
    end
    %frame=flipud(frame);
    
    %frameC = frame(y1:y2, x1:x2); % crop frames
    frameC = frame(y1:y2, :);
    %frameD = frame(:, x1:x2) %preserve the data on the row
    frameD = frame;
    %Normalizing 

    frameC64 = double(frameC) + 1; %to convert from 16 to double 
    
    original_stack(:, :, i) = double(frameD) + 1; %put the first frame on the first page
    crop_stack(:, :, i) = frameC64;
%     frameN = zeros(H,L); %matrix to store normalized figure 
% 
%     for l = 1:L 
%     
%         Max = max(frameC64(:,l)); %max value in a column 
%         Min = min(frameC64(:,l)); %min value in a column 
%         nMin = 0; %new min
%         nMax = 1; %new max 
%     
%         for h = 1:H 
%         
%              Norm = (frameC64(h,l)-Min)*((nMax-nMin)/(Max-Min))+nMin ;
%         
%             frameN(h,l)=Norm;
%         
%         end
%     
%     end

    frameN = (frameC64 - min(min(frameC64))) /(max(max(frameC64)) - min(min(frameC64)));
    
    frameN=medfilt2(frameN,[3 3]); %Applying 2D median filter% 
    frame_size = size(frameN);
    bowman_depth = zeros(1, frame_size(1, 2));
    %start = 15;%this is where to start detect bowman
    %i = 6;
    
    for k = 1:frame_size(1, 2) 
    row_array = frameN(:, k);
    value_array = smoothdata(row_array, "gaussian", 10);
    start_indice = start_offset + 5;
    top10_num = round(frame_size(1, 1) * top_ratio);%the size of top_ratio of row_array
    top10_max = maxk(row_array, top10_num);%store the top 10% number 
    avg_top10 = mean(top10_max);%mean of the average top 10%
    low40_num = round(frame_size(1, 1) * low_ratio);%the size of low_ratio of row_array
    low40_min = mink(row_array, low40_num);%store the low 40% number 
    avg_low40 = mean(low40_min);
    top2low_ratio = avg_top10 / avg_low40;
    while true
        nearby_avg_right = start_indice - start_offset;
        %nearby_avg_left = start_indice - 20;
        if start_indice > frame_size(1, 1)
            break %if reach the start of image array, then stop
        else
            if value_array(start_indice) > top2low_ratio * ...
                    mean(value_array(1 : nearby_avg_right))
                bowman_depth(1, k) = start_indice; %bowman brightness should 
                %be greater than its sourrounding area times the multiplier
                %leaving the gap between the start_indice and nearby interval
                %help reliably detect the peak
                break
            else
                start_indice = start_indice + 1; %if we don't find the brightness
                %then keep moving
            end    
        end
    end
    
    
    end


%% eliminate outlier
    isNZ=(~bowman_depth==0);           % addressing logical array of nonzero elements
    
    
    smooth_bowman = smoothdata(bowman_depth(isNZ), 'movmedian', 5);
    
    x3 = 1:frame_size(1, 2);
    P3=polyfit(x3(isNZ),smooth_bowman(1, :), 2);
    yy3_partial=polyval(P3, x3(isNZ));
    yy3_complete = polyval(P3, x3);
    
    delta=yy3_complete-bowman_depth(1, :); %computes the distance between the polynomial and the peaks
    STD_delta=std(yy3_partial - smooth_bowman);
    %avg_bowman = mean(smooth_bowman);
    for clmn=1:frame_size(1, 2)
        if (abs(delta(clmn))>STD_delta * 0.5)
            bowman_depth(1,clmn)=yy3_complete(1, clmn);%eliminate outllier by replacing it with avg value
        end
    end
    polysf(i, :) = bowman_depth;
end

%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[xq, yq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_row, cornea_column, cornea_v, xq, yq);

%% plot surface

%surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

smooth_surf = smooth2a(surface_cornea, 2, 2); %yy4 surface reconstruction  
smooth_poly = smooth2a(polysf, 15, 15);
%figure;

%surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')

%% transpose matrix
%surface_cornea = surface_cornea';
smooth_surf = smooth_surf';
%% introduce randomness in the smooth_surf to deal with artifact
%normal_random = normrnd(0, 2, size(smooth_surf));
random_smooth_surf = smooth_surf;
%% flatten the surface(based on the minimal point of the surface)
center_point = [round(poly_size(1, 1)/2), round(poly_size(1, 2)/2)];
center_surface = random_smooth_surf((center_point(1) - 200) : (center_point(1) + 200), (center_point(2) - 150) : (center_point(2) + 150));
%center_surface = center_surface;
%before find the minimum value, eliminate the edge where artifacts ususally
%occurs
peak_point = floor(min(min(center_surface)));% the level that the surface will be flattened to 
delta = round(random_smooth_surf - peak_point);%obtain the difference
delta_matrix = fillmissing(delta,'nearest');%replace possible NaN value in d
%% flatten the surface(manually choose the point that are relatively smooth)

% column_num = round(size_of_crop(2)/2);
% depth_frame = original_stack(:, column_num, :);
% depth_frame = squeeze(depth_frame);
% imshow(depth_frame, cmap);
% [~, x_cor] = ginput(1);%frame number of the selected point
% peak_point = random_smooth_surf(round(x_cor), column_num);% point of relative smooth surface
% delta = round(random_smooth_surf - peak_point);
% delta_matrix = fillmissing(delta, "nearest");
end