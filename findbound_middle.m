function [plotpeak, leftbound, rightbound] = findbound_middle(image, threshold)
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
b = eyec/2; 
leftbound = eyec/2;%start at the middle of the image
while true
    if b - 20 > 0
        if mean(plotpeak(b:(b + 20))) == 0 && plotpeak(b - 1) - plotpeak(b) > 0%check if the intensity on the right_hand side is 
            %0 and the intensity on the left_hand side is above 0
            if mean(plotpeak((b - 20):(b - 1))) > plotpeak(b - 1) / 10 %this is to check the average intensity around that point
                %if it is low then we know that this is isolated point and
                %it is not boundary
                leftbound = b;
                break
            else
                b = b - 1;
            end
        else
            b = b - 1;
        end
    else
        disp('no left bound')
        break
    end
end

%rightbound

a = eyec/2;
rightbound = eyec/2;
while true
    if a + 20 < eyec
        if mean(plotpeak((a - 20) : a)) == 0 && plotpeak(a + 1) - plotpeak(a) > 0
            if mean(plotpeak((a + 1):(a + 20))) > plotpeak(a + 1)/10
                rightbound = a;
                break
            else 
                a = a + 1;
            end
        else
            a = a + 1;
        end
    else
        disp('no right bound')
        break
    end
end


end