function [pp, dp1, leftbound, rightbound] = findbound_derivative(image)

%in the step below I normalize the intensity of the image 
Average = double(nearest(image));
normalization=(Average - min(min(Average)))/(max(max(Average)) - min(min(Average)));
[~, eyec] = size(image);
%plot
zeroplot = zeros(1, eyec);

%In the step below, I calculate the average intensity of each column in
%image normalization
for i = 1:eyec
    zeroplot(1, i) = mean(normalization(:, i));
end


pp = zeroplot;
plot = smooth(pp,0.01, 'rloess');% use 'smooth' function to smooth the curve in intensity plot
dp1 = diff(plot);% use 'diff' function to find the approximate derivative of the intensity plot
midpoint = round(length(dp1)/2);
%find the left boundary

leftdp1 = -dp1(1:midpoint);%select the left_hand side section of the image dp1 and flip the image upside down
[pk, loca] = findpeaks(leftdp1);%find the location of the local minimum intensity in dp1
[~, mi] = max(pk);
first_min_location = loca(mi);%find the position of global minimum
pk(mi) = [];
loca(mi) = [];
[~, smi] = max(pk);
second_min_location = loca(smi);%find the position of the second greatest minimum
leftbound = max([first_min_location, second_min_location]);% define the location of the left bound is the
%maximum value among the location of first and second global minimum
%intensity

%find the rightbound
rightdp1 = dp1(midpoint:(eyec - 1));%select the right_hand side setion of the image dp1
[pks, locs] = findpeaks(rightdp1);%
[~, lo] = max(pks);
first_max_location = locs(lo);
pks(lo) = [];
locs(lo) = [];
[~, la] = max(pks);
second_max_location = locs(la);
rightbound =round(midpoint + min([first_max_location, second_max_location]));%the location of right boundary is 
%the minimum of the location of first and second global maximum
%
end

    


