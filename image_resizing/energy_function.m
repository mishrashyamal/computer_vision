function energy = energy_function(img)
    gray_img = rgb2gray(img);

    sigma = 2;
    filter_size = 2*ceil(3*sigma)+1; % size of filter
    [X,Y] = meshgrid(-filter_size/2:filter_size/2, -filter_size/2:filter_size/2);
    filter_kernel = exp(-(X.^2 + Y.^2)/(2*sigma^2)) / (2*pi*sigma^2);
    filter_kernel = filter_kernel / sum(filter_kernel(:)); % normalize filter kernel

    smooth_img = conv2(gray_img, filter_kernel, 'same');

    [Gx, Gy] = imgradientxy(smooth_img);
    G = sqrt(Gx.^2 + Gy.^2);

    energy = G;
end
