function [miceye1, lower2, midd2, upper2] = cornea_layer_notnormal(path, num1, upper_depth, middle_depth, lower_depth)
%column_avg = 50 by default, column_depth = 5, middle_depth =3, lower_depth
%= 5
% path = 'G:\\OCM_MCFM_Data\\E\\Buffer_Double_Frame\\Mice_21_July_2021\\Mice1\\Eye1\\cornea5_PBS_1V\\';
%pupil1
s1 = 'image_';
% num1 = randi([1 997], 1, 1);
%  num1 = 911;
num2 = num1 + 1;
num3 = num1 + 2;
num4 = num1 + 3;
s2 = num2str(num1);%picture number 2
s2_2 = num2str(num2);
s2_3 = num2str(num3);
s2_4 = num2str(num4);
s3 = '.asc';
% path1 = strcat(path, s1, s2, s3);
% path2 = strcat(path, s1, s2_2, s3);
% path3 = strcat(path, s1, s2_2, s3);
F1 = load([path, strcat(s1, s2, s3)]);
F2 = load([path, strcat(s1, s2_2, s3)]);
F3 = load([path, strcat(s1, s2_3, s3)]);
F4 = load([path, strcat(s1, s2_4, s3)]);
eye1 = rot90(F1);
eye2 = rot90(F2);
eye3 = rot90(F3);
eye4 = rot90(F4);
Average = double(nearest(eye1));
Average2 = double(nearest(eye2));
Average3 = double(nearest(eye3));
Average4 = double(nearest(eye4));
miceye1 =(Average - min(min(Average)))/(max(max(Average)) - min(min(Average)));
miceye2 =(Average2 - min(min(Average2)))/(max(max(Average2)) - min(min(Average2)));
miceye3 =(Average3 - min(min(Average3)))/(max(max(Average3)) - min(min(Average3)));
miceye4 =(Average4 - min(min(Average4)))/(max(max(Average4)) - min(min(Average4)));
[eyer, eyec] = size(eye1);
%% upper layer
upper = zeros(2, eyec);

for i = 1:eyec
    colum_avg = (miceye1(:, i) + miceye2(:, i) + miceye3(:, i) + miceye4(:, i)) / 4;
%     colum_avg = miceye1(:, i);
%     colum_avg = (miceye1(:, i) + miceye2(:, i)) / 2;
    colum_smooth = colum_avg;
%     colum_smooth = smooth(colum_avg, 5);
    start = 16;
    while true
        if start < eyer
           if max(colum_smooth((start - 15) : (start - 5))) > 0 && upper_depth * ...
                   mean(colum_smooth((start - 15) : (start - 5)))< colum_smooth(start)
                
               new_column = [i; start];
               upper(:, i) = new_column;
               break
           else           
               start = start + 1;
           end
        else
             %disp('no upper');
             break
        end
    end
end
%% eliminate outliers
upper1 = upper(:, ~(upper(2, :) == 0));
upper_outlier_tolerant = 30;
for m = 3 : (length(upper1(1, :)) - 2)
    if (upper1(1, m) - upper1(1, m - 1))^2 + (upper1(2, m) - upper1(2, m - 1))^2 > upper_outlier_tolerant^2 && ...
             (upper1(1, m) - upper1(1, m + 1))^2 + (upper1(2, m) - upper1(2, m + 1))^2 > upper_outlier_tolerant^2
         upper1(2, m) = 0;
    end
end

% for m = 2 : (length(upper1(1, :)) - 1)
%     if (upper1(1, m) - upper1(1, m - 1))^2 + (upper1(2, m) - upper1(2, m - 1))^2 > upper_outlier_tolerant^2 || ...
%              (upper1(1, m) - upper1(1, m + 1))^2 + (upper1(2, m) - upper1(2, m + 1))^2 > upper_outlier_tolerant^2
%          upper1(2, m) = 0;
%     end
% end
% while m + 2 <= length(upper1(1, :))
%     if (upper(1, m) - upper(1, m - 1))^2 + (upper(2, m) - upper(2, m - 1))^2 > upper_outlier_tolerant^2 && ...
%             (upper(1, m) - upper(1, m + 1))^2 + (upper(2, m) - upper(2, m + 1))^2 > upper_outlier_tolerant^2
%         upper1(:, m) = [];
%     elseif m + 2 == length(upper1(1, :))
%         break
%     else
%         m = m + 1;
%     end
% end

upper2 = upper1(:, ~(upper1(2, :) == 0));
upper2(2, :) = nearest(hampel(upper2(2, :), 10));% we need to make sure that the values in upper2 are integer, since we it represents
%the index of the pixels 
%% middle layer
upper_dim = size(upper2);
dim = upper_dim(2);
midd = zeros(2, eyec);
for k = 1 : dim
    index = upper2(1, k);
    colum_avg = (miceye1(:, index) + miceye2(:, index) + miceye3(:, index) + miceye4(:, index)) / 4;
%     colum_avg = miceye1(:, i);
%     colum_avg = (miceye1(:, index) + miceye2(:, index)) / 2;
    start = upper2(2, k) + 40;
    bound = start + 80;
    while bound + 21 < eyer 
        if start <  bound 
            if middle_depth * max(colum_avg((start - 20):(start - 10))) < colum_avg(start) && mean(colum_avg((start + 10):(start + 20))) > 0.02            
                new_column = [index; start];
                midd(:, index) = new_column;
                break
            else
                start = start + 1;
            end
        else
            %disp('no middle');
            break
        end
    end
end
%% eliminate outlier
midd1 = midd(:, ~(midd(2, :) == 0));

    midd_outlier_tolerant = 30;
for m = 3 : (length(midd1(1, :)) - 2)
    if (midd1(1, m) - midd1(1, m - 2))^2 + (midd1(2, m) - midd1(2, m - 2))^2 > midd_outlier_tolerant^2 && ...
             (midd1(1, m) - midd1(1, m + 2))^2 + (midd1(2, m) - midd1(2, m + 2))^2 > midd_outlier_tolerant^2
         midd1(2, m) = 0;
    end
end
%     while k <= length(midd1)
%         if 2 * (max(midd1(2, (k - 5):(k - 1))) - min(midd1(2, (k - 5):(k - 1)))) < abs(midd1(2, k) - midd1(2, k - 1)) ...
%                 && max(midd1(2, (k - 5):(k - 1))) - min(midd1(2, (k - 5):(k - 1))) > 0%if max = min, even small variations will be
%             %deleted
%             midd1(:, k) = [];
%         elseif k == length(midd1)
%             break
%         else
%             k = k + 1;
%         end
%     end
% 
%     m = length(midd1) - 10;
%     while m >= 1 
%         if m + 10 <= length(midd1) && 2 * (max(midd1(2, (m + 1):(m + 6))) - min(midd1(2, (m + 1):(m + 6)))) < abs(midd1(2, m) - midd1(2, m + 1)) ...
%                 && max(midd1(2, (m + 1):(m + 6))) - min(midd1(2, (m + 1):(m + 6))) > 0%if max = min, even small variations will be
%             %deleted
%             midd1(:, m) = [];
%         elseif m == 1
%             break
%         else
%             m = m - 1;
%         end
%     end

midd2 = midd1(:, ~(midd1(2, :) == 0));
%% bottom layer

lower = zeros(2, eyec);

for i = 1:eyec
    colum_avg = (miceye1(:, i) + miceye2(:, i) + miceye3(:, i) + miceye4(:, i)) / 4;
%     colum_avg = miceye1(:, i);
%     colum_avg = (miceye1(:, i) + miceye2(:, i)) / 2;
    colum_smooth = flip(colum_avg);
%     colum_smooth = smooth(colum_avg, 5);
    start = 102;
    while true
        if start < eyer - 50
           if max(colum_smooth((start - 50) : (start - 10))) > 0 && lower_depth * ...
               max(colum_smooth((start - 50) : (start - 10))) < colum_smooth(start) && ...
           mean(colum_smooth((start + 5) : (start + 25))) > 0.09%first criteria of test,
           %this is to find the great peak "start" in the vertical column. the
           %great peak could potentially be the upper layer 
           
                   real_start = eyer - start;%reverse the location of start 
%                    lowlen = length(lower(1, :));
               %if lowlen > 30 && 4 * (max(lower(2, (lowlen - 30):lowlen)) - ...
                   %min(lower(2, (lowlen - 30):lowlen))) > abs(real_start - lower(2, lowlen))%second criteria of test,
               %this is to compare the location of the peak "real_start" we
               %obtain with the location of the upper layer nearby. their y
               %coordinate location should not be far away from each other
                   
                    new_column = [i; real_start];
                    lower(:, i) = new_column;
                    break
           else           
               start = start + 1;%if fail the first criteria, we proccede to the next 
           end
        else
             %disp('no bottom');
             break
        end
    end
end
%% eliminate outliers
lower1 = lower(:, ~(lower(2, :) == 0));

    lower_outlier_tolerant = 30;
for m = 2 : (length(lower1(1, :)) - 2)
    if (lower1(1, m) - lower1(1, m - 1))^2 + (lower1(2, m) - lower1(2, m - 1))^2 > lower_outlier_tolerant^2 && ...
             (lower1(1, m) - lower1(1, m + 1))^2 + (lower1(2, m) - lower1(2, m + 1))^2 > lower_outlier_tolerant^2%if the point
         %is too far away from other points nearby, then we will consider
         %it as outliers
         lower1(2, m) = 0;
    end
end

% while k <= length(lower1)
%     if 2 * (max(lower1(2, (k - 20):(k - 5))) - min(lower1(2, (k - 20):(k - 5)))) < abs(lower1(2, k) - lower1(2, k - 1)) ...
%             && max(lower1(2, (k - 20):(k - 5))) - min(lower1(2, (k - 20):(k - 5))) > 0%if max = min, even small variations will be
%         %deleted
%         lower1(:, k) = [];
%     elseif k == length(lower1)
%         break
%     else
%         k = k + 1;
%     end
% end
% 
% m = length(lower1) - 21;
% while m >= 1
%     if m + 21 <= length(lower1) && 2 * (max(lower1(2, (m + 5):(m + 20))) - min(lower1(2, (m + 5):(m + 20)))) < abs(lower1(2, m) - lower1(2, m + 1)) ...
%             && max(lower1(2, (m + 5):(m + 20))) - min(lower1(2, (m + 5):(m + 20))) > 0%if max = min, even small variations will be
%         %deleted
%         lower1(:, m) = [];
%     elseif m == 1
%         break
%     else
%         m = m - 1;
%     end
% end
lower2 = lower1(:, ~(lower1(2, :) == 0));
lower2(2, :) = hampel(lower2(2, :), 10);%use hampel function to eliminate outliers in the curve.
end
%% check

% imshow(miceye1);
% hold on
% line(lower1(1, :), lower1(2, :), 'Color', 'yellow');
% line(midd1(1, :), midd1(2, :), 'Color', 'red');
% line(upper1(1, :), upper1(2, :), 'Color', 'blue');
% hold off

%%

% cord = 100;
% colum_one = (miceye1(:, cord) + miceye2(:, cord) + miceye3(:, cord) + miceye4(:, cord)) / 4;
% colum_two = smooth(colum_one, 5);
% len = length(colum_one);
% figure;
% plot(colum_one)