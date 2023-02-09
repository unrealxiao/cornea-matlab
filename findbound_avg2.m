function [pp, plot, leftbound, rightbound] = findbound_avg2(image)
%PP is the image after average every column in the original frame, plot is
%the image after smoothing pp
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

plot = smooth(pp,0.01, 'rloess');
low_intens = min(plot);
midpoint = round(eyec/2);

%leftbound
left = midpoint - 50;%instead starting at the middle point, we start at the point from where the middle point is 50 pixels away
%the reason is that we want to make sure that the average intensity of the
%interval between 'left' and 'midpoint' is stable
while true
   if left > 0
       if plot(left) - low_intens > 3 * (mean(plot((left + 1) : midpoint)) - low_intens)%this is the test we use to determine if 
           %'left' is the boundary. the idea is that the intensity on the
           %left boundary should be much higher than the average intensity
           %on the right hand side of the left boundary.if it passes this
           %test then we believe 'left' is the boundary
           break
       else
           left = left - 1;%if it doesn't pass the test then we keep going 
       end
   else
       disp('no left bound found')
       break
   end
end
leftbound = left;

%rightbound
%the way we find the right bound is the same with the way we use to find
%the left bound
right = midpoint + 50;
while true
   if right < eyec
      if plot(right) - low_intens > 3 * (mean(plot(midpoint : (right - 1))) - low_intens)
          break
      else
          right = right + 1;
      end
   else
       disp('no right bound found')
       break
   end
end
rightbound = right;
end