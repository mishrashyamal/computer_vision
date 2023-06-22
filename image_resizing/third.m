img = imread('image10.jpg');

gray_img = rgb2gray(img);

sigma = 2;
smooth_img = imgaussfilt(gray_img, sigma);

[Gx, Gy] = imgradientxy(smooth_img);

energy_map = sqrt(Gx.^2 + Gy.^2);   

seam_mat = energy_map;

for i = 2:size(seam_mat, 1)
    for j = 1:size(seam_mat, 2)
        % Handle edge cases
        if j == 1
            seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j), seam_mat(i-1, j+1));
        elseif j == size(seam_mat, 2)
            seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j-1), seam_mat(i-1, j));
        else
            seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j-1), min(seam_mat(i-1, j), seam_mat(i-1, j+1)));
        end
    end
end

%  backtracing
[~, idx] = min(seam_mat(end, :));
opt_seam = zeros(size(seam_mat, 1), 1);
opt_seam(end) = idx;

for i = size(seam_mat, 1)-1:-1:1
    if opt_seam(i+1) == 1
        [~, m] = min([seam_mat(i, 1), seam_mat(i, 2)]);
        opt_seam(i) = m;
    elseif opt_seam(i+1) == size(seam_mat, 2)
        [~, m] = min([seam_mat(i, end-1), seam_mat(i, end)]);
        opt_seam(i) = size(seam_mat, 2) - (m-1);
    else
        [~, m] = min([seam_mat(i, opt_seam(i+1)-1), seam_mat(i, opt_seam(i+1)), seam_mat(i, opt_seam(i+1)+1)]);
        opt_seam(i) = opt_seam(i+1) + (m-2);
    end
end

img_seam = img;

for i = 1:size(opt_seam, 1)
    img_seam(i, opt_seam(i), 1) = 255;
    img_seam(i, opt_seam(i), 2) = 0;
    img_seam(i, opt_seam(i), 3) = 0;
end


figure;
subplot(1, 2, 1); imshow(img); title('Original Image');
subplot(1, 2, 2); imshow(img_seam); title('Image with Optimal Seam');
