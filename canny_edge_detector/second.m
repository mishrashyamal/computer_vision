org_img = imread('circles1.gif');

if size(org_img, 3) == 3
    org_img = rgb2gray(org_img);
end

Px = [-1 0 1; -1 0 1; -1 0 1];
Py = [-1 -1 -1; 0 0 0; 1 1 1];

K = 15; 
sigma = 2; 
h = fspecial('gaussian', K, sigma);
smoothed_I = conv2(double(org_img), double(h), 'same');

Ix = conv2(smoothed_I, Px, 'same');
Iy = conv2(smoothed_I, Py, 'same');

[abs_Ix, abs_Iy, mag_I] = deal(abs(Ix), abs(Iy), sqrt(Ix.^2 + Iy.^2));

Ix_org = conv2(double(org_img), Px, 'same');
Iy_org = conv2(double(org_img), Py, 'same');

[abs_Ix_org, abs_Iy_org, mag_I_org] = deal(abs(Ix_org), abs(Iy_org), sqrt(Ix_org.^2 + Iy_org.^2));

% Display
figure;
subplot(2, 4, 1), imshow(org_img, []), title('Original Image')
subplot(2, 4, 2), imshow(abs_Ix, []), title('Absolute value of change in x (Smoothed)')
subplot(2, 4, 3), imshow(abs_Iy, []), title('Absolute value of change in y (Smoothed)')
subplot(2, 4, 4), imshow(mag_I, []), title('Overall magnitude of the gradient (Smoothed)')
subplot(2, 4, 5), imshow(org_img, []), title('Original Image')
subplot(2, 4, 6), imshow(abs_Ix_org, []), title('Absolute value of change in x (Original)')
subplot(2, 4, 7), imshow(abs_Iy_org, []), title('Absolute value of change in y (Original)')
subplot(2, 4, 8), imshow(mag_I_org, []), title('Overall magnitude of the gradient (Original)')
