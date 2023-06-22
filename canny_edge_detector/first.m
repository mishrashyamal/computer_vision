org_img = imread('circles1.gif');

if size(org_img, 3) == 3
    org_img = rgb2gray(org_img);
end

K_values = [3, 9, 11, 15];   % 4 K values
sigma_values = [0.5, 1, 2, 4];   % 4 corresponding sigma values

figure;
subplot_idx = 1;
for k = 1:length(K_values)
    K = K_values(k);
    sigma = sigma_values(k);

    [x, y] = meshgrid(-(K-1)/2:(K-1)/2, -(K-1)/2:(K-1)/2);
    h = exp(-(x.^2 + y.^2) / (2*sigma^2));
    h = h / sum(h(:)); 

    smoothed_im = myConv(double(org_img), double(h));

    subplot(2, 4, subplot_idx); 
    imshow(org_img); 
    title('Original');
    
    subplot(2, 4, subplot_idx + 4); 
    imshow(uint8(smoothed_im)); 
    title(['K=', num2str(K), ', \sigma=', num2str(sigma)]);

    subplot_idx = subplot_idx + 1;
end

function out = myConv(org_img, kernel)
    [m, n] = size(org_img);
    [k, l] = size(kernel);
    pad_m = floor(k/2);
    pad_n = floor(l/2);
    padded_im = padarray(org_img, [pad_m, pad_n], 'replicate', 'both');
    out = zeros(m, n);
    for i = 1:m
        for j = 1:n
            out(i, j) = sum(sum(padded_im(i:i+k-1, j:j+l-1) .* kernel));
        end
    end
end
