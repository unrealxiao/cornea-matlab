function [plotpeak, leftbound, rightbound] = findbound(image, threshold, zero_distance, height)
%FINDBOUND find the boundary of pupil in an image
%image is the image we are going to process
%threshold is the variable in function peakdet.
%the bigger the gap inbetween pupil boundary, the longer zero_distance should be


Average = double(nearest(image));
normalization=(Average - min(min(Average)))/(max(max(Average)) - min(min(Average)));
[eyer, eyec] = size(image);
plotpeak = zeros(1, eyec);
%plotpeak 
for i = 1 : eyec;
    scan = normalization(:,i);
    [maxtab, mintab] = peakdet(scan, threshold);
    
    if (isempty(maxtab) == 1) %If there is NO peaks at all (show an empty matrix as a result), plot nothing.
        plotpeak(1, i) = 0;
    else
        plotpeak(1, i) = max(maxtab(:,1));%if there are peaks, plot all peaks
    end
    
 
end
%leftbound
i = 1; 
leftbound = 0;
while true;
    if i < (eyec - zero_distance);%choose 40 to make sure we obtain the correct bound
        if (plotpeak(i) - plotpeak(i + 1)) > height && mean(plotpeak((i + 1) : (i + zero_distance + 1))) == 0;
            leftbound = i;
            %leftbound;
            break
        end
        i = i + 1;
        
    else
        print('probably no left boundary');
        break
    end
end

%rightbound

a = leftbound;
rightbound = 0;
while true;
    if a < eyec;
        if (plotpeak(a + 1) - plotpeak(a)) > height && mean(plotpeak((a - zero_distance) : a)) == 0;
            rightbound = a;
            %rightbound;
            break
        end
        a = a + 1;
        
    else
        print('probably no right boundary');
        break
    end
end


end

