clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
mkdir(Path,'Flat_ENDO_part2'); % Create a subfolder in the directeory
Path_save=[Path,'Flat_ENDO_part2\'];
numframes = 1000; % number of input frames
num=1; %initial frame number 
deg=2; %degree of the polynomial fit
deg2=2;
threshold=.3; % for the peak detection 
name=[Path,name];
[I,cmap] = dicomread(name);  % find out size of images by importing one
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
name_5= 'frame';

% Perform the parallel computing

polysf=zeros(numframes,L); %yy4 surface reconstruction  

parfor i = 1:numframes

    count=i+num-1;   

    index = num2str(count);  % build filename with index
    filename = strcat(name_5,index);
    
    [frame,cmap] = dicomread([Path, filename]); % open frame
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
    P3=polyfit(x4,y4,deg);
    yy4=polyval(P3, xx2);
%     hold on; plot(xx2,yy4,'r','LineWidth',1); hold on;
%     figure
%     imshow(frameN); hold on; plot(xx2,yy4,'r','LineWidth',1); 

    polysf(i,:)=yy4;

end
    
    
%%    
%Smoothing the surface
smoothsf=smooth2a(polysf, 20, 20); %yy4 surface reconstruction  


    
%Shift 
%frameSh=zeros(H*2,L); %Matrix giving more room for the shifted regions 
const=floor(min(min(smoothsf)));
surf(smoothsf, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');

%obtain the delta matrix 

delta_matrix = round(smoothsf - const);

parfor k = 1:numframes
    
    index = num2str(k);  % build filename with index
    filename = strcat(name_5,index);
    
    [frame,~] = dicomread([Path, filename]);%read the frame
    
    F = frame(:,x1:x2); %delete the extra column
    F64 = double(F) +1 ; %convert to double
    %F64 = imgaussfilt(F64, 2);
    S2F=size(F);
    H2F=S2F(1,1);
    L2F=S2F(1,2);
    frameShtwo=zeros(size(F64)); 
    
    %moving pixel on F64
    
    for col = 1:L
        frameShtwo(1:(H2F - delta_matrix(k, col)), col) = F64((delta_matrix(k, col) + 1):H2F, col);
        frameShtwo((H2F - delta_matrix(k, col) + 1):H2F, col) = F64(1:delta_matrix(k, col), col);
    end
    
    if strcmp(Cl,'uint8')==1
        frameShtwo=uint8(frameShtwo); 
    elseif strcmp(Cl,'uint16')==1
        frameShtwo=uint16(frameShtwo); 
    end

    



    fig_index = num2str(k);  % build filename with index
    imgname = strcat(Path_save,'frameSh',index,Ext);

    dicomwrite(frameShtwo,imgname);
  
end

toc

%% eliminating artifact from enfaced image
%making one slide only

address_a = Path_save;
name_1 = 'frameSh';
name_2 = num2str(200);
name_3 = '.DCM';

[image, cmap] = dicomread([address_a, strcat(name_1, name_2, name_3)]);

[~, column] = size(image);%obtain the column size of top view image
top_view = zeros(1000, column);


depth = const; %define the depth of the image

%% stack all image together

for i = 1:1000
   name_2 = num2str(i); 
   frame = dicomread([address_a, strcat(name_1, name_2, name_3)]);
   %frame = imgaussfilt(frame, 1);
   slide = double(frame(36, :));%have to set it to double, otherwise smooth function won't work
   top_view(i, :) = slide;
end

[top_raw, top_col] = size(top_view);

