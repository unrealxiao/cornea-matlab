function [delta_matrix, smooth_surf, smooth_poly, Path_save_flat, Path_save_cross_flat, Path_save_cross_unflat, numframes, original_stack, crop_stack, flip, Cl] = cornea_delta(thres_hold, flip)
%if you want to flip the images(so that endothelium will be on top), enter
%flip as "Y"
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
%Ext=name(end-3:end);
mkdir(Path,'Flat'); % Create a subfolder in the directeory
mkdir(Path, 'cross_flat');
mkdir(Path, 'cross_unflat')
Path_save_flat=[Path,'Flat\'];
Path_save_cross_flat = [Path, 'cross_flat\'];
Path_save_cross_unflat = [Path, 'cross_unflat\'];
%count the number of files in the directories
dcm_file = dir([Path, '*.DCM']);

numframes = length(dcm_file); % number of input frames

num=1; %initial frame number 
deg=3; %degree of the polynomial fit
deg2=2;
threshold=thres_hold; % for the peak detection 
name=[Path,name];
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
x1 = min(floor(p(1)), floor(p(2))); %xmin
y1 = min(floor(p(3)), floor(p(4))); %ymin
x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
y2 = max(ceil(p(3)), ceil(p(4)));   %ymax

frameC2=frame(y1:y2, x1:x2); 
imshow(frameC2,cmap) %frameC = cropped frame 

size_of_crop = size(frameC2);

S=size(frameC2);
H=S(1,1); %height of figure (rows) 
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
    frameC = frame(y1:y2, x1:x2); % crop frames
    frameD = frame(:, x1:x2) %preserve the data on the row
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
    
    %FrameN=medfilt2(frameN,[3 3]); %Applying 2D median filter% 
    %figure
    %imshow(frameN) 

    %%Peak Detection 

    %ypixels = 1:H; % # of y-pixels = # of rows
    %frameR = zeros(H,L);

    frameR2= zeros(S);
    %xx2=1:L; 
    yy2=1:H;

    for clmn = 1:L 
    
        scan=frameN(:,clmn);
        scan=smooth(scan);                     
        [maxtab, ~]=peakdet(scan,threshold,yy2); 
   
        if (isempty(maxtab)==1) %Puts 1 if there is(are) any peak(s) detected
         [~,pos]=max(scan);  
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


    [y4, x4]=find(frameR22==1); %will be top surface of the lens
    cord_arry = zeros(1, L);
    for in = 1:length(x4)
       cord_arry(x4(in)) = y4(in); %add coordinate of the peak into the cord_arry
    end

%     [y4, x4]=find(frameR22==1); %will be top surface of the lens
%     P3=polyfit(x4,y4,deg);
%     yy4=polyval(P3, xx2);

    polysf(i,:)=cord_arry;%cord_arry contains the coordinate of the cornea surface, 0 are placed if no coordinates found
%     polysf(i,:)=yy4;
end

%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[xq, yq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_row, cornea_column, cornea_v, xq, yq);

%% plot surface

%surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

smooth_surf = smooth2a(surface_cornea, 5, 5); %yy4 surface reconstruction  
smooth_poly = smooth2a(polysf, 15, 15);
%figure;

%surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')

%% transpose matrix
%surface_cornea = surface_cornea';
smooth_surf = smooth_surf';
%% flatten the surface
center_point = [round(poly_size(1, 1)/2), round(poly_size(1, 2)/2)];
center_surface = smooth_surf((center_point(1) - 200) : (center_point(1) + 200), (center_point(2) - 150) : (center_point(2) + 150));
%center_surface = center_surface;
%before find the minimum value, eliminate the edge where artifacts ususally
%occurs
peak_point = floor(min(min(center_surface)));% the level that the surface will be flattened to 
delta = round(smooth_surf - peak_point);%obtain the difference
delta_matrix = fillmissing(delta,'nearest');%replace possible NaN value in d

end