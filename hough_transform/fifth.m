% Load the grayscale image
img = imread('circles1.gif');

if size(img, 3) == 3
    img = rgb2gray(img);
end

img = imresize(img, 0.5);

% Apply edge detection to get a binary image
binary_img = edge(img, 'canny', [0.2 0.6]);

% Apply Hough Circle detection to the binary image
rMin = 20; % lower bound of radius range
rMax = 300; % upper bound of radius range
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

% Find the centers and radii of all circles using Hough Circle detection
centers = [];
radii = [];
x0_all = [];
y0_all = [];
r_all = [];
for r = rMin:rMax
    [y,x] = find(hough_circle(:,:,r-rMin+1) >= 0.7*max(hough_circle(:)));
    centers = [centers; [x,y]];
    radii = [radii; repmat(r, length(x), 1)];
    x0_all = [x0_all; x];
    y0_all = [y0_all; y];
    r_all = [r_all; repmat(r, length(x), 1)];
end


% Remove circles that are fully enclosed in others
enclosed = false(size(radii));
for i = 1:length(radii)
    for j = i+1:length(radii)
        if sqrt(sum((centers(i,:) - centers(j,:)).^2)) + radii(i) <= radii(j)
            enclosed(i) = true;
            break;
        end
    end
end
centers(enclosed,:) = [];
radii(enclosed) = [];

% Remove circles that are partial or fully outside of the image
outside = false(size(radii));
for i = 1:length(radii)
    if centers(i,1)-radii(i) < 1 || centers(i,1)+radii(i) > size(binary_img,2) || centers(i,2)-radii(i) < 1 || centers(i,2)+radii(i) > size(binary_img,1)
        outside(i) = true;
    end
end
centers(outside,:) = [];
radii(outside) = [];



% Display the results
figure;

subplot(1,3,1); % original image
imshow(img);
title('Original Image');

subplot(1,3,2); % binary image
imshow(binary_img);
title('Binary Image');

subplot(1,3,3); % original image with major circles superimposed
imshow(img);
title('Major Circles');
viscircles(centers, radii, 'EdgeColor', 'r', 'LineWidth', 1);

disp('x0   y0   r');
for i = 1:length(x0_all)
    fprintf('%-4d %-4d %-4d\n', x0_all(i), y0_all(i), r_all(i));
end


