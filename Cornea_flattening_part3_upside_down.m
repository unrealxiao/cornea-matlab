clear variables
close all 
% Open a reference image from the stack, indetify the directory and create
% a new directory called 'Flat_ENDO' to store the results of the flattening
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
Ext=name(end-3:end);
mkdir(Path,'Flat_ENDO'); % Create a subfolder in the directeory
mkdir(Path, 'Average_Flat_ENDO_ENFACE');
Path_save=[Path,'Flat_ENDO\'];
Path_save_avg=[Path, 'Average_Flat_ENDO_ENFACE\'];
numframes = 500; % number of input frames
num=1; %initial frame number 
deg=2; %degree of the polynomial fit
deg2=2;
threshold=.3; % for the peak detection 
name=[Path,name];
[I,cmap] = dicomread(name);  % find out size of images by importing one
I=flipud(I);
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


    [y4, x4]=find(frameR22==1); %will be top surface of the lens
    cord_arry = zeros(1, L);
    for in = 1:length(x4)
       cord_arry(x4(in)) = y4(in); %add coordinate of the peak into the cord_arry
    end

%     [y4, x4]=find(frameR22==1); %will be top surface of the lens
%     P3=polyfit(x4,y4,deg);
%     yy4=polyval(P3, xx2);


%     polysf(i,:)=yy4;
    polysf(i,:)=cord_arry;
    %cord_arry contains the coordinate of the cornea surface, 0 are placed if no coordinates found
end

%% eliminate the 0 in polysf and make the surface of cornea

[cornea_row, cornea_column, cornea_v] = find(polysf); %find the value and coordinate of all non-zero value
poly_size = size(polysf);
[xq, yq] = meshgrid(1:poly_size(1, 1), 1:poly_size(1, 2));
surface_cornea = griddata(cornea_row, cornea_column, cornea_v, xq, yq, 'nearest');

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
center_surface = smooth_surf((center_point(1) - 200) : (center_point(1) + 200), (center_point(2) - 150) : (center_point(2) + 150));
%before find the minimum value, eliminate the edge where artifacts ususally
%occurs
peak_point = floor(min(min(center_surface)));% the level that the surface will be flattened to 
delta = round(smooth_surf - peak_point);%obtain the difference
delta_matrix = fillmissing(delta,'nearest');%replace possible NaN value in d

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
 %% check one cross-section
% address_a = Path_save;
% name_1 = 'frameSh';
% name_2 = num2str(800);
% name_3 = '.DCM';
% 
% [image, cmap] = dicomread([address_a, strcat(name_1, name_2, name_3)]);
% imshow(image, cmap)
% [row, column] = size(image);%obtain the column size of top view image
% %% top view
% parfor middle_frame = (row - max(y1, y2)) : (row - max(3, min(y1, y2)))
%     top_view = zeros(1000, column);
%     for i = 1:1000
%        name_2 = num2str(i); 
%        frame = dicomread([address_a, strcat(name_1, name_2, name_3)]);
%        %frame = imgaussfilt(frame, 1);
%        slide_1 = double(frame(middle_frame - 2, :));
%        slide_2 = double(frame(middle_frame - 1, :));
%        slide_3 = double(frame(middle_frame, :));
%        slide_4 = double(frame(middle_frame + 1, :));
%        slide_5 = double(frame(middle_frame + 2, :));%have to set it to double, otherwise smooth function won't work
%        top_view(i, :) = (slide_1 + slide_2 + slide_3 + slide_4 + slide_5) / 5;
%     end
% 
%     if strcmp(Cl,'uint8')==1
%         top_view=uint8(top_view); 
%     elseif strcmp(Cl,'uint16')==1
%         top_view=uint16(top_view); 
%     end
% 
%     imgname_save = strcat(Path_save_avg,'frameSh',num2str(middle_frame),'.DCM');
%     dicomwrite(top_view, imgname_save);
% end
% 
% %% top view test
% 
% middle_frame = row - (peak_point + y1 - 1);
% top_view = zeros(1000, column);
% for i = 1:1000
%    name_2 = num2str(i); 
%    frame = dicomread([address_a, strcat(name_1, name_2, name_3)]);
%    %frame = imgaussfilt(frame, 1);
%    slide_1 = double(frame(middle_frame - 2, :));
%    slide_2 = double(frame(middle_frame - 1, :));
%    slide_3 = double(frame(middle_frame, :));
%    slide_4 = double(frame(middle_frame + 1, :));
%    slide_5 = double(frame(middle_frame + 2, :));%have to set it to double, otherwise smooth function won't work
%    top_view(i, :) = (slide_1 + slide_2 + slide_3 + slide_4 + slide_5) / 5;
% end
