function img = remove_vertical_seam(img, seam)
   
    [rows, cols, ~] = size(img);
    
    output_cols = cols - 1;
    output_img = zeros(rows, output_cols, 3, class(img));
    
    for i = 1:rows
        output_img(i, :, :) = img(i, [1:seam(i)-1, seam(i)+1:end], :);
    end
    
    img = output_img;
end
