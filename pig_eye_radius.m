%close all; clear all;
%%set up picture
path = 'E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\pupildata\\test_data\\prediction\\';
%pupil1
s1 = 'prediction_';
s2 = '309';%picture number 2
s2_2 = '310';%picture number 3
s3 = '.png';
path1 = strcat(path, s1, s2, s3);
path2 = strcat(path, s1, s2_2, s3);
F1 = load([path, strcat(s1, s2, s3)]);
F2 = load([path, strcat(s1, s2_2, s3)]);
eye1 = rot90(F1);
eye2 = rot90(F2);
Average = double(nearest(eye1));
Average2 = double(nearest(eye2));
pigeye1 =(Average - min(min(Average)))/(max(max(Average)) - min(min(Average)));
pigeye2 = (Average2 - min(min(Average2)))/(max(max(Average2)) - min(min(Average2)));
[eyer, eyec] = size(eye1);
%% find the plotpeak of pigeye1
[p1, plt, left, right] = findbound_avg2(eye1);

[eyer, eyec] = size(eye1);
figure;
imshow(pigeye1);
x = left;
y = 0:eyer;
hold on
plot(x*ones(size(y)), y)

x2 = right;
hold on
plot(x2*ones(size(y)), y)

% scatter(1:eyec, p1)
% hold off
%the first image looks fine
%% find the pupil boundary distance in pigeye1
% tic
% [~, pic, left2, right2] = findbound_avg2(eye1);
% toc
% [eyer, eyec] = size(eye1);
plot(plt);
hold on
xline(x)
hold on
xline(x2)


hold off

% p2 = smooth(p2,0.1, 'rloess');
% figure;
% plot(pic)
% hold on
% xline(x)
% hold on
% xline(x2)

%scatter(1:eyec, p2);

% midpoint = round(length(p2)/2);
% rightdp2 = p2(midpoint:length(p2));
% [pks, locs] = findpeaks(rightdp2);
% [~, lo] = max(pks);
% first_max_location = locs(lo);
%hold off

%adjust the threshhold from 0.6 to 0.4, and then we got the correct result

%%
%findbound_avg
tic
[pp, p2, left2, right2] = findbound_derivative(eye1);
toc
[eyer, eyec] = size(pigeye1);
imshow(pigeye1);
x = left2;
y = 0:eyer;
hold on
plot(x*ones(size(y)), y)

x2 = right2;
hold on
plot(x2*ones(size(y)), y);


hold off

% p2 = smooth(p2,0.1, 'rloess');
figure;
plot(pp)
hold on
xline(x)
hold on
xline(x2)
figure;
plot(p2)

hold on
xline(x)
hold on
xline(x2)


hold off

% figure;


% plot(dp1)
% [pks, locs] = findpeaks(dp1(round(954/2):953));
% [p, l] = max(pks);
% pks(l) = [];
% locs(l) = [];
% [pa, la] = max(pks);

%%
%for eye1
sset = zeros(4, 1000);
path = 'E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\Pig_Eyeball-2021-6-29\\Eye1\\pupil2\\';
s1 = 'image_';
s3 = '.asc';
parfor i = 1:1000
    [radius, n1, n2] = findbound_random(path, s1, s3);
    dat = zeros(4, 1);
    dat(1, 1) = i; %row 1 stores the index
    dat(2, 1) = radius; % row 2 stores the radius
    dat(3, 1) = n1; % row 3 stores the frame number
    dat(4, 1) = n2; 
    sset(:, i) = dat;
    i
end
%%
scatter(sset(3, :), sset(2, :))
xlim([200 700])
ylim([2 3.5])
title('radius vs frame index');
xlabel('frame index');
ylabel('radius');
%%

