clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
numframes = 1000; % number of input frames
num=1; %initial frame number 
deg=2; %degree of the polynomial fit
deg2=2;
threshold=.3; % for the peak detection 
name=[Path,name];
[I,cmap] = dicomread(name);  % find out size of images by importing one
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
name_5= 'frame';

% Perform the parallel computing

polysf=zeros(numframes,L);

parfor i = 1:numframes

    count=i;
    %count = 300;
    %index = sprintf('%04d' ,count);
    index = num2str(count);% build filename with index
    filename = strcat(name_5,index, '.DCM');
    
    [frame,~] = dicomread([Path, filename]);
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

    %FrameN=medfilt2(frameN,[3 3]); %Applying 2D median filter% 
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
    [y34, x34]=find(frameR22==1); %will be top surface of the lens
    cord_arry = zeros(1, L);
    for in = 1:length(x34)
       cord_arry(x34(in)) = y34(in) + min(y1, y2); %add coordinate of the peak into the cord_arry
    end
    polysf(i,:)=cord_arry;
end
 %% detect bottom
 
 frame = I;
imshow(frame,cmap) 
p = ginput(2);  % have user crop image by selecting 2 coordinate points 
p(p<1)=1;
%x1 = min(floor(p(1)), floor(p(2))); %xmin
yb = min(floor(p(3)), floor(p(4))); %ymin
%x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
yb2 = max(ceil(p(3)), ceil(p(4)));   %ymax

frameC2=frame(yb:yb2, x1:x2); 
imshow(frameC2,cmap) %frameC = cropped frame 

S=size(frameC2);
H=S(1,1); %height of figure (rows) 
L=S(1,2); %lenght of figure (columns)


name_5= 'frame';

% Perform the parallel computing

polysf2=zeros(numframes,L); %yy4 surface reconstruction  

parfor i = 1:numframes

    count=i;
    %count = 300;
    %index = sprintf('%04d' ,count);
    index = num2str(count);% build filename with index
    filename = strcat(name_5,index, '.DCM');
    
    [frame,~] = dicomread([Path, filename]);
    %frame = rot90(frame,0); % rotate frames, this can be used if the frames are not in the proper orientation
    %frame=flipud(frame);
    frameC = frame(yb:yb2, x1:x2); % crop frames
    
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

    %FrameN=medfilt2(frameN,[3 3]); %Applying 2D median filter% 
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
    [y4, x4]=find(frameR2==1); %will be top surface of the lens
    P4=polyfit(x4,y4,deg);
    yy4=polyval(P4, x4);
%     hold on; plot(x3,yy3,'b','LineWidth',1); hold on;
%     figure
%     imshow(frameN); hold on; plot(xx2,yy3,'r','LineWidth',1);


%% Refine the polynomial fit by removing bad peaks

    delta=yy4-y4; %computes the distance between the polynomial and the peaks
    STD_delta=std(delta);
    frameR22=frameR2;
    
    for clmn=1:L
        if (abs(delta(clmn))>STD_delta/0.9)
            frameR22(:,clmn)=0;
        end
    end


%     imshow(frameR22) 
    [y44, x44]=find(frameR22==1); %will be top surface of the lens
    cord_arry = zeros(1, L);
    for in = 1:length(x44)
       cord_arry(x44(in)) = y44(in) + min(yb, yb2); %add coordinate of the peak into the cord_arry
    end
%%end
%  figure
%  imshow(frame, cmap)
%  hold on
%  line(x34, y34 + min(y1, y2), 'Color', 'red');
%  line(x44, y44 + min(yb, yb2), 'Color', 'yellow');
 polysf2(i,:)=cord_arry;%cord_arry contains the coordinate of the cornea surface, 0 are placed if no coordinates found


end

%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[xq, yq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_row, cornea_column, cornea_v, xq, yq);

[cornea_row2, cornea_column2, cornea_v2] = find(polysf2); %find the value and coordinate of all non-zero value
poly_size2 = size(polysf2);
[xq2, yq2] = meshgrid(1:poly_size2(1, 1), 1:poly_size2(1, 2));
surface_cornea2 = griddata(cornea_row2, cornea_column2, cornea_v2, xq2, yq2);

%% smooth and plot the surface

smooth_surf = smooth2a(surface_cornea, 5, 5); %yy4 surface reconstruction  
smooth_surf2 = smooth2a(surface_cornea2, 20, 20);

%%
figure;

surf(smooth_surf * 0.7, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
hold on
surf(smooth_surf2 * 0.7, 'FaceColor','y', 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'ZDir','reverse')
set(gca, 'xtick', [], 'ytick', [])
hold off
%% depth map

epi_stroma = smooth_surf2 - smooth_surf;
figure;
imagesc(epi_stroma * 0.7);% convert to micro meter
% caxis([0 80]);
title('epi');
figure;