%% 
close all; clear all;
eye = imread('/Users/mac/Desktop/lab/pupil/eye.png', 'png');
eye2 = eye(280:321, 94:755);%segmentation. the section of pupil

Average = double(nearest(eye2));
normalization=(Average - min(min(Average)))/(max(max(Average)) - min(min(Average)));
[eyer, eyec] = size(eye2);
plotpeak = zeros(1, eyec);
for i = 1 : eyec;
    scan = normalization(:,i);
    [maxtab, mintab] = peakdet(scan, 0.6);
    
    if (isempty(maxtab) == 1) %If there is NO peaks at all (show an empty matrix as a result), plot nothing.
        plotpeak(1, i) = 0;
    else
        plotpeak(1, i) = max(maxtab(:,1));%if there are peaks, plot all peaks
    end
    
 
end

x = 1 : eyec;
scatter(x, plotpeak);

%next we are going to find the boundary of the pupil 
%find the left boundary of the pupil
%left_candidate = zeros(1, eyec);

%% left boundary
i = 1; 
leftbound = 0;
while true;
    if i < (eyec - 40);
        if (plotpeak(i) - plotpeak(i + 1)) > 15 && mean(plotpeak((i + 1) : (i + 41))) == 0;
            leftbound = i;
            leftbound
            break
        end
        i = i + 1;
        
    else
        print('probably no left boundary');
        break
    end
end

%then we need to find the right boundary if the pupil

%% right boundary
a = leftbound;
rightbound = 0;
while true;
    if a < (eyec - 40);
        if (plotpeak(a + 1) - plotpeak(a)) > 15 && mean(plotpeak((a - 30) : a)) == 0;
            rightbound = a;
            rightbound
            break
        end
        a = a + 1;
        
    else
        print('probably no right boundary');
        break
    end
end

%% check the boundary line on the image

imshow(eye2)
x = leftbound;
y = 0:eyer;
hold on
plot(x*ones(size(y)), y)

x2 = rightbound;
hold on
plot(x2*ones(size(y)), y)

%% stretch the image eye2
[rows, columns, nCC] = size(eye2);
%subplot(2, 1, 1);
%imshow(eye2);
newcolumns = round(1.5 * columns);
%subplot(2, 1, 2);
stretchedImage = imresize(eye2, [rows newcolumns]);
imshow(stretchedImage);
figure, imshow(eye2);

%% 

Average_1 = double(nearest(stretchedImage));
normalization=(Average_1 - min(min(Average_1)))/(max(max(Average_1)) - min(min(Average_1)));
[eyer_1, eyec_1] = size(stretchedImage);
plotpeak_1 = zeros(1, eyec_1);
for i = 1 : eyec_1;
    scan = normalization(:,i);
    [maxtab, mintab] = peakdet(scan, 0.6);
    
    if (isempty(maxtab) == 1) %If there is NO peaks at all (show an empty matrix as a result), plot nothing.
        plotpeak_1(1, i) = 0;
    else
        plotpeak_1(1, i) = max(maxtab(:,1));%if there are peaks, plot all peaks
    end
    
 
end

x = 1 : eyec_1;
scatter(x, plotpeak_1);

%% left bound
i = 1; 
leftbound_1 = 0;
while true;
    if i < (eyec - 40);
        if (plotpeak_1(i) - plotpeak_1(i + 1)) > 15 && mean(plotpeak_1((i + 1) : (i + 41))) == 0;
            leftbound_1 = i;
            leftbound_1
            break
        end
        i = i + 1;
        
    else
        print('probably no left boundary');
        break
    end
end

%% right bound

b = leftbound_1;
rightbound_1 = 0;
while true;
    if b < (eyec_1 - 40);
        if (plotpeak_1(b + 1) - plotpeak_1(b)) > 15 && mean(plotpeak_1((b - 30) : b)) == 0;
            rightbound_1 = b;
            rightbound_1
            break
        end
        b = b + 1;
        
    else
        print('probably no right boundary');
        break
    end
end

%% check

imshow(stretchedImage)
x = leftbound_1;
y = 0:eyer_1;
hold on
plot(x*ones(size(y)), y)

x2 = rightbound_1;
hold on
plot(x2*ones(size(y)), y)