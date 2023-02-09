function[radius, num1, num2] = findbound_random(path, s1, s3)
%FINDBOUND_RANDOM randomly select two image from folder and compute the
%radius
num1 = randi([200 700], 1, 1);
num2 = num1 + 1;
s2 = num2str(num1);
s2_2 = num2str(num2);
% path1 = strcat(path, s1, s2, s3);
% path2 = strcat(path, s1, s2_2, s3);
% pigeye1 = imread(path1);
% pigeye2 = imread(path2);
F1 = load([path, strcat(s1, s2, s3)]);
F2 = load([path, strcat(s1, s2_2, s3)]);
pigeye1 = rot90(F1);
pigeye2 = rot90(F2);
% [~, ~, left1, right1] = findbound_derivative(pigeye1);
% [~, ~, left2, right2] = findbound_derivative(pigeye2);
[~, ~, left1, right1] = findbound_avg2(pigeye1);
[~, ~, left2, right2] = findbound_avg2(pigeye2);
% [~, left1, right1] = findbound_middle(pigeye1, 0.4);
% [~, left2, right2] = findbound_middle(pigeye2, 0.4);
% pigeye1(:, left1) = 60000;
% pigeye1(:, right1) = 60000;
% pigeye2(:, left2) = 60000;
% pigeye2(:, right2) = 60000;
% image1 = imfuse(pigeye1, pigeye2, 'montage');
gap1 = (right1 - left1)*0.0108;
gap2 = (right2 - left2)*0.0108;
radius = diam_2(gap1, gap2, 1.62);









