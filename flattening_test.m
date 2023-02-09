clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.tif'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
mkdir(Path,'Flat_ENDO_part2'); % Create a subfolder in the directeory
Path_save=[Path,'Flat_ENDO_part2\'];
numframes = 1000; % number of input frames
num=1; %initial frame number  %degree of the polynomial fit
deg2=2;
threshold=.25; % for the peak detection 
name=[Path,name];
[I,cmap] = imread(name);  % find out size of images by importing one
I=flipud(I);
Cl=class(I);
tic

%% 

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
deg=3;
name_5= 'frame';

% Perform the parallel computing

%     polysf=zeros(numframes,L); %yy4 surface reconstruction  
% 
%     i = 500;
% 
%     count=i;   
% 
%     index = num2str(count);  % build filename with index
%     filename = strcat(name_5,index);
%     
%     
%     
%     [frame,cmap] = dicomread([Path, filename]); % open frame%remember to
% %     uncomment this 
% 
%  
%     
%     
%     
%     frame=flipud(frame);
    frame = I;
    frameC = frame(y1:y2, x1:x2); % crop frames
    
    %Normalizing 

    frameC64 = double(frameC); %to convert from 16 to double 
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

    FrameN = medfilt2(frameN,[3 3]); %Applying 2D median filter% 
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
   
        if (isempty(maxtab)==1); %Puts 1 if there is(are) any peak(s) detected
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

%% Refine the polynomial fit by removing bad peaks from endo flattening upside down

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
    imshow(frameN); hold on; plot(xx2,yy4,'r','LineWidth',1); 

%% Refine the polynomial fit from cornea_flattening_part 3 upside down

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
    imshow(frameN); hold on; plot(x4,y4,'r','LineWidth',1); 

    
 %%  flatten one frame only
    
    F = frame(:,x1:x2); %delete the extra column
    %[F, cmap] = imread(name);
    %F = F(:, x1:x2);
    F64 = double(F) +1 ; %convert to double
    F64 = flipud(F64);%flip the images
    %F64 = imgaussfilt(F64, 0.8);
    S2F=size(F);
    H2F=S2F(1,1);
    L2F=S2F(1,2);
    frameShtwo=zeros(size(F64)); 
    
    %moving pixel on F64
    delta_matrix = round(max(yy4) - yy4);
    k = 1;
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
    %frameShtwo = flipud(frameShtwo);
    