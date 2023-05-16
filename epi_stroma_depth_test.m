clear
close all
%% check which frame has clear stroma boundary
Path=[pwd,'/'];
[name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 

dcm_file = dir([Path, '*.DCM']);
%%
[I,cmap] = dicomread([Path, 'frame160.dcm']);

imshow(I,cmap) 

p = ginput(2);  % have user crop image by selecting 2 coordinate points 
p(p<1)=1;
x1 = min(floor(p(1)), floor(p(2))); %xmin
y1 = min(floor(p(3)), floor(p(4))); %ymin
x2 = max(ceil(p(1)), ceil(p(2)));   %xmax
y2 = max(ceil(p(3)), ceil(p(4)));   %ymax
%%
[I,cmap] = dicomread([Path, 'frame250.dcm']);
frameC2=I(y1:y2, :);
frameC2 = double(frameC2) + 1;
frameC2 = (frameC2 - min(min(frameC2))) /(max(max(frameC2)) - min(min(frameC2)));
frameC2 = imadjust(frameC2);
imshow(frameC2, cmap)

%% find the boundary of bowman layer on every column
start_offset = 20;%start counting from the middle, instead of from the end to avoid 
%noise on the boundary of the images
bowman_to_nearby_ratio = 1.6;%ratio of bowman brightness to nearby surrounding
close all
frame_med = medfilt2(frameC2);
frame_size = size(frame_med);
bowman_depth = zeros(1, frame_size(1, 2));
start = frame_size(1, 1) - start_offset;%this is where to start detect bowman
if start < 0
   disp('start lower than 0 for the first place. need to be > 0')
   return
end
%i = 6;

for i = 1:frame_size(1, 2) 
row_array = frame_med(:, i);
value_array = smoothdata(row_array, "gaussian", 10);
start_indice = start;

while true
    nearby_avg_left = start_indice + 10;
    nearby_avg_right = start_indice + 20;
    if start_indice <= 0
        break %if reach the start of image array, then stop
    elseif nearby_avg_left > frame_size(1, 1)
        disp(['nearby_avg_left exceeds row bound, error at i', num2str(i)])
        return % can't exceed the length of image array
    elseif nearby_avg_right > frame_size(1, 1)
        disp(['nearby_avg_right exceeds row bound, error at i', num2str(i)])
        return % can't exceed the length of image array
    else
        if value_array(start_indice) > bowman_to_nearby_ratio * ...
                mean(value_array(nearby_avg_left : nearby_avg_right))
            bowman_depth(1, i) = start_indice; %bowman brightness should 
            %be greater than its sourrounding area
            break
        else
            start_indice = start_indice - 1; %if we don't find the brightness
            %then keep moving
        end    
    end
end


end
plot(value_array)
xline(start_indice, "color", "red");
figure
imshow(frame_med,cmap)
%% eliminate outlier
x3 = 1:frame_size(1, 2);
P3=polyfit(x3,bowman_depth(1, :), 2);
yy3=polyval(P3, x3);

delta=yy3-bowman_depth(1, :); %computes the distance between the polynomial and the peaks
STD_delta=std(delta);
avg_bowman = mean(bowman_depth);
for clmn=1:frame_size(1, 2)
    if (abs(delta(clmn))>STD_delta * 0.7)
        bowman_depth(1,clmn)=avg_bowman;%eliminate outllier by replacing it with avg value
    end
end

isNZ=(~bowman_depth==0);           % addressing logical array of nonzero elements

figure
imshow(frame_med,cmap)
hold on
smooth_curve = smoothdata(bowman_depth(isNZ), 'movmedian', 5);
plot(flipud(smooth_curve), 'color', 'yellow')
%plot(flipud(bowman_depth))