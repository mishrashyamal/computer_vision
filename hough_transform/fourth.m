% Load the grayscale image
img = imread('circles1.gif');


% Apply edge detection to get a binary image
binary_img = edge(img, 'canny', [0.2 0.6]);

% Apply Hough Circle detection to the binary image
rMin = floor(min(size(binary_img))/6); % lower bound of radius range
rMax = floor(max(size(binary_img))/3); % upper bound of radius range


numRadii = rMax - rMin + 1;
hough_circle = zeros(size(binary_img, 1), size(binary_img, 2), numRadii);

for r = rMin:rMax
    for y0 = 1:size(binary_img,1)
        for x0 = 1:size(binary_img,2)
            if binary_img(y0,x0) == 1 
                for theta = 0:pi/180:2*pi 
                    x = round(x0 + r*cos(theta)); 
                    y = round(y0 + r*sin(theta)); 
                    if x > 0 && x <= size(binary_img,2) && y > 0 && y <= size(binary_img,1)
                        hough_circle(y,x,r-rMin+1) = hough_circle(y,x,r-rMin+1) + 1;
                    end
                end
            end
        end
    end
end

% Find the dominant circle
[maxVal, maxIdx] = max(hough_circle(:));
[y0Idx, x0Idx, rIdx] = ind2sub(size(hough_circle), maxIdx);
x0 = x0Idx;
y0 = y0Idx;
r = rIdx + rMin - 1;

% Display the results
figure;

subplot(1,3,1); % original image
imshow(img);
title('Original Image');

subplot(1,3,2); % binary image
imshow(binary_img);
title('Binary Image');

subplot(1,3,3); % original image with dominant circle superimposed
imshow(img);
title(sprintf('Original image with dominant circle superimposed in red.'));
viscircles([x0,y0], r, 'EdgeColor', 'r', 'LineWidth', .5);
