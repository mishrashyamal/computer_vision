image = imread('test-image.jpg');

grayscale = 0.2989*image(:,:,1) + 0.5870*image(:,:,2) + 0.1140*image(:,:,3);

norm_img = im2double(grayscale);

binary1 = norm_img >= 0.25;
binary2 = norm_img >= 0.5;
binary3 = norm_img >= 0.75;

figure;
subplot(2,3,1), imshow(image), title('Color Image')
subplot(2,3,2), imshow(grayscale), title('Grayscale Image')
subplot(2,3,3), imshow(binary1), title('t=25%')
subplot(2,3,4), imshow(binary2), title('t=50%')
subplot(2,3,5), imshow(binary3), title('t=75%')
