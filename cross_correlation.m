function [total_shift] = cross_correlation(stack)
size_stack = size(stack);

%correlation_correct_stack = zeros(size_stack);%creat another empty stack for correlation frame.

number_of_frames = size_stack(1, 3);

shift_indice = zeros(1, number_of_frames - 1); %store the shift for each frames

%we will determine the shift between all adjacent frames, and then store
%the shift in the shift_indice array. the shift_indice(1) will store the
%shift between flatten_stack(:, :, 1) and flatten_stack(:, :, 2)
stack_2 = stack; %create another flatten stack to avoid overhead in parloop
parfor j = 1 : (number_of_frames - 1)
    FF1 = stack(:, :, j);
    FF2 = stack_2(:, :, j + 1);
    l1=sum(FF1);
    %l2=sum(FF2);
    len=length(l1);
    
    %corel=xcorr(l1,l2);

    %figure; plot (corel);grid
    %[M, shift]=max(corel)

    Mcorel=xcorr2(FF1,FF2);
    [~, pos]=max(max(Mcorel))
    shift_indice(1, j)=pos-len;
end

%compute the total shift for each frame
total_shift = zeros(1, number_of_frames - 1);

for n = 1 : (number_of_frames - 1)
    total_shift(1, n) = sum(shift_indice(1, 1 : n));%total_shift store the total shift for each frames
end
end