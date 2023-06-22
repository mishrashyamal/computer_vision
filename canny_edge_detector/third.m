image = imread('circles1.gif');

if size(image, 3) == 3
    image = rgb2gray(image);
end

kernel_size = 9;
sigma = 1;

gaussian_kernel = fspecial('gaussian', kernel_size, sigma);

smoothed_image = conv2(double(image), double(gaussian_kernel), 'same');

[Gx, Gy] = gradient(smoothed_image);

gradient_magnitude = sqrt(Gx.^2 + Gy.^2);
gradient_orientation = atan2(Gy, Gx);

gradient_orientation(gradient_orientation < 0) = gradient_orientation(gradient_orientation < 0) + 2*pi;

suppressed = zeros(size(gradient_magnitude));
for i = 2:size(gradient_magnitude, 1)-1
    for j = 2:size(gradient_magnitude, 2)-1
        theta = gradient_orientation(i,j);
        if (theta < -7*pi/8) || (theta >= 7*pi/8) || (-pi/8 <= theta && theta < pi/8)
            r1 = i;
            c1 = j+1;
            r2 = i;
            c2 = j-1;
        elseif (theta >= -5*pi/8 && theta < -3*pi/8) || (theta >= 3*pi/8 && theta < 5*pi/8)
            r1 = i-1;
            c1 = j+1;
            r2 = i+1;
            c2 = j-1;
        elseif (theta >= -3*pi/8 && theta < -pi/8) || (theta >= 5*pi/8 && theta < 7*pi/8)
            r1 = i-1;
            c1 = j;
            r2 = i+1;
            c2 = j;
        elseif (theta >= -7*pi/8 && theta < -5*pi/8) || (theta >= 1*pi/8 && theta < 3*pi/8)
            r1 = i-1;
            c1 = j-1;
            r2 = i+1;
            c2 = j+1;
        end
        
        if (gradient_magnitude(i,j) > gradient_magnitude(r1,c1) && gradient_magnitude(i,j) > gradient_magnitude(r2,c2))
            suppressed(i,j) = gradient_magnitude(i,j);
        end
    end
end

figure;
subplot(2,2,1); imshow(image); title('Original Image');
subplot(2,2,2); imshow(suppressed, []); title('Non-Maximum Suppression');
