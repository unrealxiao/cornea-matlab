%% test cornea_layer
path = 'G:\\OCM_MCFM_Data\\E\\Buffer_Double_Frame\\4th_Aug2021_mice\\Mice_Female_5_TC_Smad4FF_PBS_1V\\';
% num1 = 1;
% [miceye1, lower2, midd2, upper2] = cornea_layer(path, num1, 5, 3, 5);
 lower_3D = zeros(1000, 1000);
 midd_3D = zeros(1000, 1000);
 upper_3D = zeros(1000, 1000);
%% lower_3D aquire
parfor i = 10:997
    [~, lower2, midd2, upper2] = cornea_layer_normal(path, i, 4, 2, 4);
    lower_array = zeros(1, 1000);
    midd_array = zeros(1, 1000);
    upper_array = zeros(1, 1000);
    k = 1;
    while k <= length(lower2(1, :))
        if length(lower2(1, :)) < 10
            break
        else
            lower_array(1, lower2(1, k)) = lower2(2, k);
            k = k + 1;
        end
    end
    j = 1;
    while j <= length(midd2(1, :))
       if length(midd2(1, :)) < 10%if the size of the midd is too small, it doesn't worth the time to plot the data
           break
       else
          midd_array(1, midd2(1, j)) = midd2(2, j);
          j = j + 1;
       end    
    end
    h = 1;
    while h <= length(upper2(1, :))
       if length(upper2(1, :)) < 10
           break
       else
           upper_array(1, upper2(1, h)) = upper2(2, h);
           h = h + 1;
       end
    end
    lower_3D(i, :) = lower_array;
    midd_3D(i, :) = midd_array;
    upper_3D(i, :) = upper_array;
end
%% only run this part once
[lower_row, lower_column, lower_v] = find(lower_3D);
[midd_row, midd_column, midd_v] = find(midd_3D);%zero indicatesno value, we should ignore all zeros value in the 2D data
[upper_row, upper_column, upper_v] = find(upper_3D);%extract all nonzero value and their indices in the matrix
%% flip the data(only run this once)
lower_v = abs(600 - lower_v);
midd_v = abs(600 - midd_v);
upper_v = abs(600 - upper_v);

%% fit the data by using griddata
[xq, yq] = meshgrid(1:1000, 1:1000);
lower_vq = griddata(lower_row, lower_column, lower_v, xq, yq);
midd_vq = griddata(midd_row, midd_column, midd_v, xq, yq);
upper_vq = griddata(upper_row, upper_column, upper_v, xq, yq);
%% eliminate all the data on the edge of the 2D data_matrix(so that the artifacts around the edge can be eliminated)

for i = 1:1000
    for k = 1:1000
       if (i - 500)^2 + (k - 500)^2 > 300^2
           lower_vq(i, k) = NaN;
           midd_vq(i, k) = NaN;
           upper_vq(i, k) = NaN;
       end
    end
end

%% use smooth function

smooth_lower = smooth2a(lower_vq, 20, 20);
smooth_midd = smooth2a(midd_vq, 20, 20);
smooth_upper = smooth2a(upper_vq, 20, 20);

%% depth map
epi_stroma = smooth_upper - smooth_midd;
epi_endo = smooth_upper - smooth_lower;
stroma_endo = smooth_midd - smooth_lower;

%% plot the data

% mesh(epi_stroma, 'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','none')

% mesh(smooth_lower, 'FaceAlpha', 0.5)
% hold on
% mesh(smooth_midd, 'FaceAlpha', 0.5)
% mesh(smooth_upper, 'FaceAlpha', 0.5)
% hold off

surf(smooth_lower, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');
hold on
surf(smooth_midd, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
surf(smooth_upper,'FaceColor','b', 'FaceAlpha',0.5, 'EdgeColor','none')
hold off

%% plot the depth
imagesc(epi_stroma * 0.7);caxis([0 110]);
title('epi');
figure;
imagesc(stroma_endo * 0.7);caxis([0 110]);
title('stroma');
figure;
imagesc(epi_endo * 0.7);caxis([0 110]);
title('whole');

%% mean and std

real_epi_stroma = epi_stroma * 0.7;
real_epi_endo = epi_endo * 0.7;
real_stroma_endo = stroma_endo * 0.7;
%epi
thickness_epi = mean(real_epi_stroma, 'all', 'omitnan');%it seems like using mean(mean(A)) will give slightly different 
%result compared with mean(A, 'all')
std_epi = std(real_epi_stroma, 0, 'all', 'omitnan');
%stroma
thickness_stroma = mean(real_stroma_endo, 'all', 'omitnan');
std_stroma = std(real_stroma_endo, 0, 'all', 'omitnan');
%whole cornea
thickness_whole = mean(real_epi_endo, 'all', 'omitnan');
std_whole = std(real_epi_endo, 0, 'all', 'omitnan');