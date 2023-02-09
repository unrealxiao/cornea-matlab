%% for neural network only
%diam = xlsread('E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\Pig_Eyeball-2021-6-8\\Eye1\output\\adap_th_diameter');
diam = xlsread('E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\Xiao\\pupillary_Xiao2');
diam_radius = (diam(:, 4) - diam(:, 3));%convert the dimension back to original image size;
%diam_radius = diam(:, 3);
radius_set = zeros(2, 999);

%%
for i = 1:999
   gap1 = diam_radius(i)*0.02;
   gap2 = diam_radius(i + 1)*0.02;
   result = zeros(2, 1);
   result(1, 1) = i;
   result(2, 1) = diam_2(gap1, gap2, 1.62);
   radius_set(:, i) = result;%compute the diameter
   i;
end

%% graph

scatter(radius_set(1, :), radius_set(2, :))
%scatter(diam(:, 2), diam_radius)
% xlim([200 700])
ylim([2 5])
title('pupillary_xiao1');
xlabel('frame index');
ylabel('radius');

%% for jung_sung code

diam = xlsread('E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\Xiao\\graph_result_from_jong_sung\\diameter');
diam_radius = diam(:, 3);
radius_set = zeros(2, 999);
for i = 1:999
   gap1 = diam_radius(i)*0.02;
   gap2 = diam_radius(i + 1)*0.02;
   result = zeros(2, 1);
   result(1, 1) = i;
   result(2, 1) = diam_2(gap1, gap2, 1.62);
   radius_set(:, i) = result;% compute the diameter
   i;
end

scatter(radius_set(1, :), radius_set(2, :))
% xlim([200 700])
 ylim([2 5])
title('Thetan');
xlabel('frame index');
ylabel('radius');

%% avg2

sset = zeros(4, 1000);
path = 'E:\\3Dimage\\Intensity\\OCT_Pupillary_Test\\pupildata\\test_data\\station_pupil\\';
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
ylim([2 3])
title('radius vs frame index');
xlabel('frame index');
ylabel('radius');