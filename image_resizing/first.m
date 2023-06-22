img = imread('image10.jpg');
[h, w, ~] = size(img);

h_new = input('Enter the desired height: ');
w_new = input('Enter the desired width: ');

% Nearest Neighbor Interpolation
target_img_nn = uint8(zeros(h_new, w_new, 3));

for i = 1:h_new
    for j = 1:w_new
        x_prime = j;
        y_prime = i;
        x = round(x_prime * (w / w_new));
        y = round(y_prime * (h / h_new));

        if x < 1
            x = 1;
        elseif x > w
            x = w;
        end
        if y < 1
            y = 1;
        elseif y > h
            y = h;
        end

        target_img_nn(i, j, :) = img(y, x, :);
    end
end

% Bilinear Interpolation
new_img_bi = uint8(zeros(h_new, w_new, 3));

scale_x = (w-1) / (w_new-1);
scale_y = (h-1) / (h_new-1);

for y_prime = 1:h_new
    for x_prime = 1:w_new
       
        x = (x_prime-1) * scale_x + 1;
        y = (y_prime-1) * scale_y + 1;
        
        x1 = floor(x);
        x2 = ceil(x);
        y1 = floor(y);
        y2 = ceil(y);
        
        if x2 > w
            x2 = w;
        end
        if y2 > h
            y2 = h;
        end
        
        pixel1 = double(img(y1,x1,:));
        pixel2 = double(img(y1,x2,:));
        pixel3 = double(img(y2,x1,:));
        pixel4 = double(img(y2,x2,:));
        
        x_frac = x - x1;
        y_frac = y - y1;
        
        interpolated_pixel = (1-x_frac)*(1-y_frac)*pixel1 + x_frac*(1-y_frac)*pixel2 ...
            + (1-x_frac)*y_frac*pixel3 + x_frac*y_frac*pixel4;
        
        new_img_bi(y_prime,x_prime,:) = uint8(interpolated_pixel);
    end
end

figure;
subplot(1, 3, 1); imshow(img); title('Original Image');
subplot(1, 3, 2); imshow(target_img_nn); title('Nearest Neighbor Interpolation');
subplot(1, 3, 3); imshow(new_img_bi); title('Bilinear Interpolation');