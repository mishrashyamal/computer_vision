image = imread('test-image.jpg');

image = im2double(image);

red = image(:,:,1);
green = image(:,:,2);
blue = image(:,:,3);

gamma1 = cat(3, red.^0.2, green.^0.2, blue.^0.2);
gamma2 = cat(3, red.^1, green.^1, blue.^1);
gamma3 = cat(3, red.^50, green.^50, blue.^50);

subplot(1,4,1), imshow(image), title('Orginal Image')
subplot(1,4,2), imshow(gamma1), title('Gamma 0.2')
subplot(1,4,3), imshow(gamma2), title('Gamma 1')
subplot(1,4,4), imshow(gamma3), title('Gamma 50')
