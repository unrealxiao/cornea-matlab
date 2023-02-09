clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
if ~exist([Path, 'Flat_ENDO'], 'dir')
    mkdir(Path, 'Flat_ENDO')
end
% Create a subfolder in the directeory
%mkdir(Path, 'Average_Flat_ENDO_ENFACE');
Path_save=[Path,'Flat_ENDO\'];
%Path_save_avg=[Path, 'Average_Flat_ENDO_ENFACE\'];
numframes = 1000; % number of input frames
num=1; %initial frame number 
deg=2; %degree of the polynomial fit
deg2=2;
threshold=.2; % for the peak detection 
name=[Path,name];
[I,cmap] = dicomread(name);  % find out size of images by importing one
I=flipud(I);
Cl=class(I);
tic

disp('tic done')

%% crop

frame = I;

imshow(frame,cmap) 
p = ginput(2);  % have user crop image by selecting 2 coordinate points 
p(p<1)=1;
x1 = min(floor(p(1)), floor(p(2))); %xmin
y1 = min(floor(p(3)), floor(p(4))); %ymin
x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
y2 = max(ceil(p(3)), ceil(p(4)));   %ymax

frameC2=frame(y1:y2, x1:x2); 
imshow(frameC2,cmap) %frameC = cropped frame 

S=size(frameC2);
H=S(1,1); %height of figure (rows) 
L=S(1,2); %lenght of figure (columns) 

disp('L=S(1,2) done')

%%
name_5= 'frame';

% Perform the parallel computing

polysf=zeros(numframes,L); %yy4 surface reconstruction  

parfor i = 1:numframes

    count=i+num-1;
    %count = 300;

    index = num2str(count);  % build filename with index
    filename = strcat(name_5,index);
    
    [frame,~] = dicomread([Path, filename]); % open frame
    %frame = rot90(frame,0); % rotate frames, this can be used if the frames are not in the proper orientation
    %frame=flipud(frame);
    frame=flipud(frame);
    frameC = frame(y1:y2, x1:x2); % crop frames
    
    %Normalizing 

    frameC64 = double(frameC) + 1; %to convert from 16 to double 
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

    FrameN=medfilt2(frameN,[3 3]); %Applying 2D median filter% 
    %figure
    %imshow(frameN) 

    %%Peak Detection 

    ypixels = 1:H; % # of y-pixels = # of rows
    frameR = zeros(H,L);

    frameR2= zeros(S);
    xx2=1:L; 
    yy2=1:H;

    for clmn = 1:L 
    
        scan=frameN(:,clmn);
        scan=smooth(scan);                     
        [maxtab, mintab]=peakdet(scan,threshold,yy2); 
   
        if (isempty(maxtab)==1) %Puts 1 if there is(are) any peak(s) detected
         [Mvalue,pos]=max(scan);  
         frameR2(pos,clmn)=1; 
        else
        frameR2(maxtab(1,1),clmn)=1;
        end
    
    end

%     figure
%     imshow(frameR2) 
    [y3, x3]=find(frameR2==1); %will be top surface of the lens
    P3=polyfit(x3,y3,deg);
    yy3=polyval(P3, x3);
%     hold on; plot(x3,yy3,'b','LineWidth',1); hold on;
%     figure
%     imshow(frameN); hold on; plot(xx2,yy3,'r','LineWidth',1);

%% Refine the polynomial fit by removing bad peaks

    delta=yy3-y3; %computes the distance between the polynomial and the peaks
    STD_delta=std(delta);
    frameR22=frameR2;
    
    for clmn=1:L
        if (abs(delta(clmn))>STD_delta*0.9)
            frameR22(:,clmn)=0;
        end
    end

%     figure
%     imshow(frameR22) 
    [y4, x4]=find(frameR22==1); %will be top surface of the lens
    cord_arry = zeros(1, L);
    for in = 1:length(x4)
       cord_arry(x4(in)) = y4(in); %add coordinate of the peak into the cord_arry
    end
%%end

% imshow(frameN)
% hold on
% line(x4, y4, 'Color', 'yellow');
    polysf(i,:)=cord_arry;%cord_arry contains the coordinate of the cornea surface, 0 are placed if no coordinates found
end

disp('polysf done')

%polysf : first row correspond to the curve in the first B-scan frame
%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[yq, xq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_column, cornea_row, cornea_v, xq, yq);

disp('surface_cornea done')
%% plot surface

surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

smooth_surf = smooth2a(surface_cornea, 25, 25); %yy4 surface reconstruction  

figure;

surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% transpose matrix
surface_cornea = surface_cornea';
smooth_surf = smooth_surf';
% flatten the surface
center_point = [round(poly_size(1, 1)/2), round(poly_size(1, 2)/2)];
center_surface = smooth_surf((center_point(1) - 200) : (center_point(1) + 200), (center_point(2) - 200) : (center_point(2) + 200));
%before find the minimum value, eliminate the edge where artifacts ususally
%occurs
peak_point = floor(min(min(center_surface)));% the level that the surface will be flattened to 
delta = round(smooth_surf - peak_point);%obtain the difference
delta_matrix = fillmissing(delta,'nearest');%replace possible NaN value in d

disp('delta done')

%%
parfor k = 1:numframes
    index = num2str(k);  % build filename with index
    filename = strcat(name_5,index);
    
    [frame,~] = dicomread([Path, filename]);%read the frame
    
    F = frame(:,x1:x2); %delete the extra column
    F64 = double(F) +1 ; %convert to double
    F64 = flipud(F64);%flip the images
    %F64 = imgaussfilt(F64, 0.8);
    S2F=size(F);
    H2F=S2F(1,1);
    L2F=S2F(1,2);
    frameShtwo=zeros(size(F64)); 
    
    %moving pixel on F64
        
    for col = 1:L
        if delta_matrix(k, col) > 0
            frameShtwo(1:(H2F - delta_matrix(k, col)), col) = F64((delta_matrix(k, col) + 1):H2F, col);
            frameShtwo((H2F - delta_matrix(k, col) + 1):H2F, col) = F64(1:delta_matrix(k, col), col);
        else
            positive_delta = abs(delta_matrix(k, col));
            frameShtwo(1:positive_delta, col) = F64((H2F - positive_delta + 1):H2F, col);
            frameShtwo((positive_delta + 1):H2F, col) = F64(1:(H2F - positive_delta), col);
        end 
    end
    
  
    
    
    
    if strcmp(Cl,'uint8')==1
        frameShtwo=uint8(frameShtwo); 
    elseif strcmp(Cl,'uint16')==1
        frameShtwo=uint16(frameShtwo); 
    end

    

    %frameShtwo = imgaussfilt(frameShtwo, 0.7);%use gaussian filter to smooth the small discontinuities caused by flattening.
    %frameShtwo = medfilt2(frameShtwo, [2, 2]);%using moving average
    frameShtwo = flipud(frameShtwo);%flip the images back 

    fig_index = num2str(k);  % build filename with index
    imgname = strcat(Path_save,'frameSh',index,Ext);

    dicomwrite(frameShtwo,imgname);
  
end

disp('OCM done')

%% green FM image

disp('start green chanel')
Path=[pwd,'/'];%choose the folder of the green channel
[name,Path]=uigetfile([Path,'*.png'],' choose green channel folder and select a random image'); 
Ext=name(end-3:end);
if ~exist([Path, 'Flat_FM_green'], 'dir')
    mkdir(Path, 'Flat_FM_green')
end

%mkdir(Path, 'Average_Flat_ENDO_ENFACE');
Path_save=[Path,'Flat_FM_green\'];
image_list = dir([Path, '*.png']);
number_of_images = length(image_list);
[sample_image, cmap] = imread([Path, name]);
C2 = class(sample_image);
[sample_r, ~] = size(sample_image);

FM_delta = round(delta_matrix / 4); %account for the height scaling used in FM

disp('FM_delta done')
%% green FM stacking to 3D images

% FM_3D_block = zeros(sample_r, L, number_of_images); %column dimension are dimension of the corped images
% %sample_r should be 1000. it's the number of B-scan
% parfor n = 1 : number_of_images
%    image_name = [num2str(n), 'Green_Channel.png'];
%    image = imread([Path, image_name]);
%    image = image(:,x1:x2);
%    FM_3D_block(:, :, n) = image;
% end
% 
% disp('FM_3D_block done')
%% green FM flattening and saving 
pixel_m = 0; % when you want to move all pixels up or down by certain amount, you can use this parameter, set it to 0 when only let FM_delta move pixels
parfor k = 1:sample_r
    
%     Bscan = zeros(number_of_images, L);
%     
%    for n = 1 : number_of_images
%        %disp(num2str(n))
%        image_name = [num2str(n), 'Green_Channel.png'];
%        image = imread([Path, image_name]);
%        image = image(k,x1:x2);
%        Bscan(n , :) = image;
%    end

    Bscan = squeeze(FM_3D_block(k, :, :))';% use transpose matrix to correct the dimension after reshape
    Bscan = flipud(Bscan);
    
    
    FM_Bscan = zeros(size(Bscan));

 
    %I don't know if the orientation of the FM_delta is the same as the
    %FM_Bscan. The dimension is correct, but the orientation might not. it
    %seems like the images are correctely flattened, but if the flattening
    %is not correct, try flip the FM_Bscan to get the correct orientation
    for col = 1:L
        if FM_delta(k, col) > 0
            FM_Bscan(1:(number_of_images - FM_delta(k, col)), col) = Bscan((FM_delta(k, col) + 1):number_of_images, col);
            FM_Bscan((number_of_images - FM_delta(k, col) + 1):number_of_images, col) = Bscan(1:FM_delta(k, col), col);
        else
            positive_delta = abs(FM_delta(k, col));
            FM_Bscan(1:positive_delta, col) = Bscan((number_of_images - positive_delta + 1):number_of_images, col);
            FM_Bscan((positive_delta + 1):number_of_images, col) = Bscan(1:(number_of_images - positive_delta), col);
        end 
    end
    
  
    
    
    
    if strcmp(C2,'uint8')==1
        FM_Bscan=uint8(FM_Bscan); 
    elseif strcmp(C2,'uint16')==1
        FM_Bscan=uint16(FM_Bscan); 
    end

    


    FM_Bscan = flipud(FM_Bscan);%flip the images back 
    
    also_Bscan = FM_Bscan;
    
    %shift the pixels further as needed
    
    also_Bscan(1:(number_of_images - pixel_m), :) = FM_Bscan((pixel_m + 1):number_of_images, :);
    also_Bscan((number_of_images - pixel_m + 1):number_of_images, :) = FM_Bscan(1:pixel_m, :);
    
    
    also_Bscan = ind2rgb(also_Bscan, cmap);%convert the images from grayscale to RGB, so that the image looks like green

    fig_index = num2str(k);  % build filename with index
    imgname = strcat(Path_save,'FM',num2str(k), Ext);

    imwrite(also_Bscan,imgname);
  
end

disp('green chanel done')

%% red channel

disp('start red channel')
Path=[pwd,'/'];%choose the folder of the green channel
[name,Path]=uigetfile([Path,'*.png'],' choose red channel folder and select a random image'); 
Ext=name(end-3:end);
if ~exist([Path, 'Flat_FM_red'], 'dir')
    mkdir(Path, 'Flat_FM_red')
end

%mkdir(Path, 'Average_Flat_ENDO_ENFACE');
Path_save=[Path,'Flat_FM_red\'];
image_list = dir([Path, '*.png']);
number_of_images = length(image_list);
[sample_image, cmap] = imread([Path, name]);
C2 = class(sample_image);
[sample_r, ~] = size(sample_image);

FM_3D_block = zeros(sample_r, L, number_of_images); %column dimension are dimension of the corped images
%sample_r should be 1000. it's the number of B-scan
parfor n = 1 : number_of_images
   image_name = ['Red_Channel', num2str(n), '.png'];%change the name
   image = imread([Path, image_name]);
   image = image(:,x1:x2);
   FM_3D_block(:, :, n) = image;
end

disp('FM_3D_block done')

%% red image flatten
 % use this when you want to move the pixels in the images further more.
parfor  k = 1 : sample_r
    Bscan = squeeze(FM_3D_block(k, :, :))';% use transpose matrix to correct the dimension after reshape
    Bscan = flipud(Bscan);
    FM_Bscan = zeros(size(Bscan));

 
    %I don't know if the orientation of the FM_delta is the same as the
    %FM_Bscan. The dimension is correct, but the orientation might not. it
    %seems like the images are correctely flattened, but if the flattening
    %is not correct, try flip the FM_Bscan to get the correct orientation
    for col = 1:L
        if FM_delta(k, col) > 0
            FM_Bscan(1:(number_of_images - FM_delta(k, col)), col) = Bscan((FM_delta(k, col) + 1):number_of_images, col);
            FM_Bscan((number_of_images - FM_delta(k, col) + 1):number_of_images, col) = Bscan(1:FM_delta(k, col), col);
        else
            positive_delta = abs(FM_delta(k, col));
            FM_Bscan(1:positive_delta, col) = Bscan((number_of_images - positive_delta + 1):number_of_images, col);
            FM_Bscan((positive_delta + 1):number_of_images, col) = Bscan(1:(number_of_images - positive_delta), col);
        end 
    end
    
  
    
    
    
    if strcmp(C2,'uint8')==1
        FM_Bscan=uint8(FM_Bscan); 
    elseif strcmp(C2,'uint16')==1
        FM_Bscan=uint16(FM_Bscan); 
    end

    

    %frameShtwo = imgaussfilt(frameShtwo, 0.7);%use gaussian filter to smooth the small discontinuities caused by flattening.
    %frameShtwo = medfilt2(frameShtwo, [2, 2]);%using moving average
    FM_Bscan = flipud(FM_Bscan);%flip the images back 
    
    also_Bscan = zeros(size(FM_Bscan));
    
    %shift the pixels further as needed
    
    also_Bscan(1:pixel_m, :) = FM_Bscan((number_of_images - pixel_m + 1):number_of_images, :);
    also_Bscan((pixel_m + 1):number_of_images, :) = FM_Bscan(1:(number_of_images - pixel_m), :);
    
    %
    
    also_Bscan = ind2rgb(also_Bscan, cmap);%convert the images from grayscale to RGB, so that the image looks like green

    fig_index = num2str(k);  % build filename with index
    imgname = strcat(Path_save,'FM',num2str(k), Ext);

    imwrite(also_Bscan,imgname);
  
end

disp('red chanel done')