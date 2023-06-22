
img = imread('test-image.jpg');

img = im2double(img);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

gray = 0.2989*R + 0.5870*G + 0.1140*B;

imshow(gray);