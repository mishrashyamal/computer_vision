img = imread('test-image.jpg');

img = im2double(img);


R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
maxVal = max(max(R, G), B);
minVal = min(min(R, G), B);
delta = maxVal - minVal;
 hue = zeros(size(maxVal));
idx = (maxVal == R);
 hue(idx) = mod((G(idx) - B(idx)) ./ delta(idx), 6);
idx = (maxVal == G);
 hue(idx) = mod(((B(idx) - R(idx)) ./ delta(idx) + 2), 6); 
idx = (maxVal == B);
 hue(idx) = mod(((R(idx) - G(idx)) ./ delta(idx) + 4), 6);

 hue =  hue / 6 * 360;
S = zeros(size(maxVal));
S(maxVal > 0) = delta(maxVal > 0) ./ maxVal(maxVal > 0);
V = maxVal;

 hue =  hue + 50;
 hue( hue > 360) =  hue( hue > 360) - 360;
 hue( hue < 0) =  hue( hue < 0) + 360;


C = V .* S;
X = C .* (1 - abs(mod( hue / 60, 2) - 1));
m = V - C;
R = zeros(size(maxVal));
G = zeros(size(maxVal));
B = zeros(size(maxVal));
idx = (0 <=  hue) & ( hue < 60);
R(idx) = C(idx);
G(idx) = X(idx);
idx = (60 <=  hue) & ( hue < 120);
R(idx) = X(idx);
G(idx) = C(idx);
idx = (120 <=  hue) & ( hue < 180);
G(idx) = C(idx);
B(idx) = X(idx);
idx = (180 <=  hue) & ( hue < 240);
G(idx) = X(idx);
B(idx) = C(idx);
idx = (240 <=  hue) & ( hue < 300);
R(idx) = X(idx);
B(idx) = C(idx);
idx = (300 <=  hue) & ( hue < 360);
R(idx) = C(idx);
B(idx) = X(idx);
R = R + m;
G = G + m;
B = B + m;


img_hue_shifted = cat(3, R, G, B);
figure;
subplot(2, 2, 1); imshow(img); title('Original Image');
subplot(2, 2, 2); imshow(img_hue_shifted); title('Adjusted shifted');