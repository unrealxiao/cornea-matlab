function [bowman_depth] = cornea_layer_notnormal(path, num1, y1, y2, start_offset, bowman_to_nearby_ratio)

%start_offset, start counting from the middle, instead of from the end to avoid 
%noise on the boundary of the images
%bowman_to_nearby_ratio. ratio of bowman brightness to nearby surrounding.
%this is used to help identify the bowman. If potential layer / neaby brightness >
%this ratio then this layer is viewd as bowman
s1 = 'frame';
s2 = num2str(num1);%picture number 2
s3 = '.DCM';
[I,~] = dicomread([path, s1, s2, s3]);
frameC2=I(y1:y2,:);
frameC2 = double(frameC2) + 1;
frameC2 = (frameC2 - min(min(frameC2))) /(max(max(frameC2)) - min(min(frameC2)));
%% find the boundary of bowman layer on every column
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
            %be greater than its sourrounding area times the multiplier
            %leaving the gap between the start_indice and nearby interval
            %help reliably detect the peak
            break
        else
            start_indice = start_indice - 1; %if we don't find the brightness
            %then keep moving
        end    
    end
end


end
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


end 