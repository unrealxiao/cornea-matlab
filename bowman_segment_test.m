clear
close all
%% check which frame has clear stroma boundary
% Path=[pwd,'/'];
% [name,Path]=uigetfile([Path,'*.DCM'],' Choose a reference image in the format DCM. '); 
Path = ['F:\OK3_V1_2023-4-11\scan15\pic', '\'];

% dcm_file = dir([Path, '*.DCM']);

[I,cmap] = dicomread([Path, 'frame350.dcm']);
%%
[height, width] = size(I);
y1 = height - 230;
y2 = height - 30;

%%
% [I,cmap] = dicomread([Path, 'frame250.dcm']);
frameC2=I(y1:y2, :);
frameC2 = double(frameC2) + 1;
frameC2 = (frameC2 - min(min(frameC2))) /(max(max(frameC2)) - min(min(frameC2)));
frameC2 = imadjust(frameC2);
imshow(frameC2, cmap)

%% find the boundary of bowman layer on every column
start_offset = 20;%start counting from the middle, instead of from the end to avoid 
%noise on the boundary of the images
bowman_to_nearby_ratio = 1.4;%ratio of bowman brightness to nearby surrounding
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
%% eliminate outlier by using histcounts
%use histcounts to count the number of occurence at different depth
%from previous step. The depth that has the highest occurence should be
%the actual bowman depth location. The other depth is false detection
[N, edges] = histcounts(bowman_depth, 15);
N_len = length(N);
%[N_sort, N_index] = sort(N, 'descend');%sort the number of occurence in descending order
[~, Max_Index_N] = max(N);
%now compute the ratio of number of values that fall into the most populus
%bin, if not above 70 percent, include the adjacent two bins and
%repeat the process until the raito of number of values exceed 70 percent

Num_of_bin_up = 0;%start the additional number of bin with 0
Num_of_bin_down = 0;%start the additional number of bin with 0
bowman_depth_len = length(bowman_depth);%total number of data points
while true
    Included_occurence = N(Max_Index_N - Num_of_bin_down:Max_Index_N + Num_of_bin_up);%number of occurence included
    ratio_of_included = sum(Included_occurence) / bowman_depth_len;%compute the ratio
    if ratio_of_included < 0.7 %when the number ratio is small, it's likely that we don't have enought bowman points
        if Max_Index_N + Num_of_bin_up + 1 <= N_len && N(Max_Index_N + Num_of_bin_up + 1) / bowman_depth_len > 0.1
            Num_of_bin_up = Num_of_bin_up + 1;%if there is additional bin on the right hand side and the number of points inside
            %that bin is large enought, then it's likely that it contains
            %the true bowman depth
        elseif Max_Index_N - Num_of_bin_down - 1 >= 1 && N(Max_Index_N - Num_of_bin_down - 1) / bowman_depth_len > 0.1
            Num_of_bin_down = Num_of_bin_down + 1;%the same logic applies here
        else
            break
        end
    else
        break%if the included bin has more than 70% of the total numbers, then we stop 
        %since it's likely that we included the bins that contain all the
        %bowman depth
    end
end


upper_limit = edges(Max_Index_N + Num_of_bin_up + 1); %upper limit of bin interval
lower_limit = edges(Max_Index_N - Num_of_bin_down); %lower limit of bin interval


%now we've obtained the interval that convers the real depth of bowman
%layer, we will only retain data points that lie inside the interval. then
%use interpolation to replace the deleted data points.

bowman_coordinate = zeros(2, bowman_depth_len);%contain coordinate of bowman
bowman_coordinate(1, :) = 1:bowman_depth_len;%first row contains the x axis
bowman_coordinate(2, :) = bowman_depth;%second row contains the y axis
bowman_filter = bowman_depth >= lower_limit & bowman_depth <= upper_limit;
bowman_noise_free = bowman_coordinate(:, bowman_filter);%eliminat most noise
%now the obvious spikes have all been removed. we will use polynomial curve
%to further remove the small spikes.
P3=polyfit(bowman_noise_free(1, :), bowman_noise_free(2, :), 2);
yy3=polyval(P3, bowman_noise_free(1, :));
delta=abs(yy3-bowman_noise_free(2, :)); %computes the distance between the polynomial and the peaks
STD_delta=std(delta);
spike_filter = delta >= STD_delta * 0.7;%detect additional spikes that fall out of std range

bowman_no_spike = bowman_noise_free(:, spike_filter);%retain the data points that are inside the std range


bowman_noise_free_interp = interp1(bowman_no_spike(1, :), ...
                                   bowman_no_spike(2, :), ...
                                   bowman_coordinate(1, :));%interpolate the deleted data point



complete_curve = fillmissing(bowman_noise_free_interp, 'nearest');%fill the
%missing data NaN in the interpolated array



% isNZ=(~nwe_bowman==0); % addressing logical array of nonzero elements

figure
imshow(frame_med,cmap)
hold on
smooth_curve = smoothdata(complete_curve, 'movmedian', 10);

plot(smooth_curve, 'color', 'yellow')
% plot(bowman_noise_free(1, :), bowman_noise_free(2, :), 'Color', 'yellow');
% plot(bowman_depth, 'Color', 'red')
hold off
% figure;
% histogram(bowman_depth, 15)
%plot(flipud(bowman_depth))