function [bowman_depth] = cornea_layer_notnormal(path, num1, y1, y2, lower_sd, upper_sd, row_frac)
%only detect bowman layer, upper_sd and lower_sd define the range of the
%brightness that we are looking for. row_frac define the row vector's
%indice that are believed to store the value of bowman layer depth
s1 = 'frame';
s2 = num2str(num1);%picture number 2
s3 = '.DCM';
[I,~] = dicomread([path, s1, s2, s3]);
frameC2=I(y1:y2,:);
frameC2 = double(frameC2) + 1;
frameC2 = (frameC2 - min(min(frameC2))) /(max(max(frameC2)) - min(min(frameC2)));
%% find the boundary of bowman layer on every column
frame_size = size(frameC2);
bowman_depth = zeros(1, frame_size(1, 2));
for i = 1:frame_size(1, 2)
    value_array = frameC2(:, i);
    Avg = mean(value_array);
    sd = std(value_array);

    [row, ~] = find(value_array > Avg + lower_sd*sd & value_array < Avg + upper_sd*sd);
    row_size = size(row);
    if row_size(1, 1) < 3 && i == 1
       bowman_depth(1, i) = 0; %when you can't find anything, put 0 in it
    elseif row_size(1, 1) < 3
        bowman_depth(1, i) = bowman_depth(1, i - 1); %when you don't have a lot, copy previous data
    else
        half_len = round(row_frac * row_size(1, 1));
        bowman_depth(1, i) = row(half_len, 1);% use something in the middle
    end

end
%% eliminate outlier
x3 = 1:frame_size(1, 2);
P3=polyfit(x3,bowman_depth(1, :), 2);
yy3=polyval(P3, x3);

delta=yy3-bowman_depth(1, :); %computes the distance between the polynomial and the peaks
STD_delta=std(delta);

for clmn=1:frame_size(1, 2)
    if (abs(delta(clmn))>STD_delta * 0.9)
        bowman_depth(1,clmn)=0;
    end
end

%% add the croped rows back

bowman_depth = bowman_depth;
end 