img = imread('test-image.jpg');

gray_img = rgb2gray(img);

gray_hist = zeros(1, 256);
red_hist = zeros(1, 256);
green_hist = zeros(1, 256);
blue_hist = zeros(1, 256);

for i = 1:size(img, 1)
    for j = 1:size(img, 2)
        gray_hist(gray_img(i,j) + 1) = gray_hist(gray_img(i,j) + 1) + 1;
        red_hist(img(i,j,1) + 1) = red_hist(img(i,j,1) + 1) + 1;
        green_hist(img(i,j,2) + 1) = green_hist(img(i,j,2) + 1) + 1;
        blue_hist(img(i,j,3) + 1) = blue_hist(img(i,j,3) + 1) + 1;
    end
end

figure;
subplot(2,2,1);
bar(gray_hist);
title('Grayscale Histogram');
subplot(2,2,2);
bar(red_hist);
title('Red Channel Histogram');
subplot(2,2,3);
bar(green_hist);
title('Green Channel Histogram');
subplot(2,2,4);
bar(blue_hist);
title('Blue Channel Histogram');
