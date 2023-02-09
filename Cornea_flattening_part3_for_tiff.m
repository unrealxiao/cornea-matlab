clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.tif'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
mkdir(Path,'Flat_ENDO_part3'); % Create a subfolder in the directeory
Path_save=[Path,'Flat_ENDO_part3\'];
numframes = 1000; % number of input frames
num=1; %initial frame number 
deg=2; %degree of the polynomial fit
deg2=2;
threshold=.3; % for the peak detection 
name=[Path,name];
[I,cmap] = imread(name);  % find out size of images by importing one
%I=flipud(I);
Cl=class(I);
tic

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

%%
name_5= 'median';

% Perform the parallel computing

polysf=zeros(numframes,L); %yy4 surface reconstruction  

parfor i = 1:numframes

    count=i;
    %count = 300;
    index = sprintf('%04d' ,count);  % build filename with index
    filename = strcat(name_5,index, '.tif');
    
    [frame,~] = imread([Path, filename]); % open frame
    %frame = rot90(frame,0); % rotate frames, this can be used if the frames are not in the proper orientation
    %frame=flipud(frame);
    frameC = frame(y1:y2, x1:x2); % crop frames
    
    %Normalizing 

    frameC64 = double(frameC) + 1; %to convert from 16 to double 
    frameN = zeros(H,L); %matrix to store normalized figure 

    for l = 1:L 
    
        Max = max(frameC64(:,l)); %max value in a column 
        Min = min(frameC64(:,l)); %min value in a column 
        nMin = 0; %new min
        nMax = 1; %new max 
    
        for h = 1:H 
        
             Norm = (frameC64(h,l)-Min)*((nMax-nMin)/(Max-Min))+nMin ;
        
            frameN(h,l)=Norm;
        
        end
    
    end

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
        if (abs(delta(clmn))>STD_delta/0.9)
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

%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[xq, yq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_row, cornea_column, cornea_v, xq, yq);

%% plot surface

surf(surface_cornea, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

smooth_surf = smooth2a(surface_cornea, 5, 5); %yy4 surface reconstruction  

figure;

surf(smooth_surf, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
%% transpose matrix
surface_cornea = surface_cornea';
smooth_surf = smooth_surf';
%% flatten the surface
center_point = [round(poly_size(1, 1)/2), round(poly_size(1, 2)/2)];
center_surface = surface_cornea((center_point(1) - 200) : (center_point(1) + 200), (center_point(2) - 200) : (center_point(2) + 200));
%before find the minimum value, eliminate the edge where artifacts ususally
%occurs
peak_point = floor(min(min(center_surface)));% the level that the surface will be flattened to 
delta = round(surface_cornea - peak_point);%obtain the difference
delta_matrix = fillmissing(delta,'nearest');%replace possible NaN value in d

%%
parfor k = 1:numframes
    index = sprintf('%04d', k);  % build filename with index
    filename = strcat(name_5,index, '.tif');
    
    [frame,~] = imread([Path, filename]);%read the frame
    
    F = frame(:,x1:x2); %delete the extra column
    F64 = double(F) +1 ; %convert to double
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

    fig_index = num2str(k);  % build filename with index
    imgname = strcat(Path_save,'frameSh',num2str(k),'.DCM');

    dicomwrite(frameShtwo,imgname);
  
end
 %% check one cross-section
address_a = Path_save;
name_1 = 'frameSh';
name_2 = num2str(200);
name_3 = '.DCM';

[image, cmap] = dicomread([address_a, strcat(name_1, name_2, name_3)]);

[~, column] = size(image);%obtain the column size of top view image
%% top view
parfor i = 1:1000
   name_2 = num2str(i); 
   frame = dicomread([address_a, strcat(name_1, name_2, name_3)]);
   %frame = imgaussfilt(frame, 1);
   slide = double(frame(113, :));%have to set it to double, otherwise smooth function won't work
   top_view(i, :) = slide;
end

[top_raw, top_col] = size(top_view);