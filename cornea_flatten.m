function[flatten_stack] = cornea_flatten(unflatten_stack, delta_matrix, flip, Cl)
size_stack = size(unflatten_stack);
row_stack = size_stack(1, 1);
column_stack = size_stack(1, 2);
flatten_stack = zeros(size_stack);%create another empty stack for flattened frame
parfor k = 1:size_stack(1, 3)
    
    original_frame = unflatten_stack(:, :, k);
    frameShtwo=zeros(row_stack, column_stack); 
    
    %moving pixel on original_frame
        
    for col = 1:column_stack
        if delta_matrix(k, col) > 0
            frameShtwo(1:(row_stack - delta_matrix(k, col)), col) = original_frame((delta_matrix(k, col) + 1):row_stack, col);
            frameShtwo((row_stack - delta_matrix(k, col) + 1):row_stack, col) = original_frame(1:delta_matrix(k, col), col);
        else
            positive_delta = abs(delta_matrix(k, col));
            frameShtwo(1:positive_delta, col) = original_frame((row_stack - positive_delta + 1):row_stack, col);
            frameShtwo((positive_delta + 1):row_stack, col) = original_frame(1:(row_stack - positive_delta), col);
        end 
    end
    
  
    
    
    
    if strcmp(Cl,'uint8')==1
        frameShtwo=uint8(frameShtwo); 
    elseif strcmp(Cl,'uint16')==1
        frameShtwo=uint16(frameShtwo); 
    end

    

    %frameShtwo = imgaussfilt(frameShtwo, 0.7);%use gaussian filter to smooth the small discontinuities caused by flattening.
    %frameShtwo = medfilt2(frameShtwo, [2, 2]);%using moving average
    if (flip == "Y")
       frameShtwo = flipud(frameShtwo);
    end
    
    flatten_stack(:, :, k) = frameShtwo;%put flattened frame into flattened stack

  
end
%%random nonsense
end