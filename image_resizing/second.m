img = imread('image10.jpg');
gray_img = rgb2gray(img);

% Create a 2D Gaussian filter kernel
sigma = 2;
filter_size = 2*ceil(3*sigma)+1; % size of filter
[X,Y] = meshgrid(-filter_size/2:filter_size/2, -filter_size/2:filter_size/2);
filter_kernel = exp(-(X.^2 + Y.^2)/(2*sigma^2)) / (2*pi*sigma^2);
filter_kernel = filter_kernel / sum(filter_kernel(:));

smooth_img = conv2(gray_img, filter_kernel, 'same');

%  gradient magnitude 
[Gx, Gy] = imgradientxy(smooth_img, 'sobel');
energy_map = sqrt(Gx.^2 + Gy.^2);

figure;
subplot(1, 2, 1); imshow(img); title('Original Image');
subplot(1, 2, 2); imshow(energy_map, []); title('Energy Map');

fprintf('Gaussian filter sigma: %d\n', sigma);
