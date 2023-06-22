
im = imread('circles1.gif');

if size(im, 3) == 3
    im = rgb2gray(im);
end

K = 9;
sigma = 1;

h = fspecial('gaussian', K, sigma);

smoothed_im = conv2(double(im), double(h), 'same');

[Gx, Gy] = gradient(smoothed_im);

mag = sqrt(Gx.^2 + Gy.^2);
ang = atan2(Gy, Gx);

ang(ang < 0) = ang(ang < 0) + 2*pi;

suppressed = zeros(size(mag));
for i = 2:size(mag, 1)-1
    for j = 2:size(mag, 2)-1
        theta = ang(i,j);
        if (theta < -7*pi/8) || (theta >= 7*pi/8) || (-pi/8 <= theta && theta < pi/8)
            r1 = i;
            c1 = j+1;
            r2 = i;
            c2 = j-1;
        elseif (theta >= -5*pi/8 && theta < -3*pi/8) || (theta >= 3*pi/8 && theta < 5*pi/8)
            r1 = i-1;
            c1 = j+1;
            r2 = i+1;
            c2 = j-1;
        elseif (theta >= -3*pi/8 && theta < -pi/8) || (theta >= 5*pi/8 && theta < 7*pi/8)
            r1 = i-1;
            c1 = j;
            r2 = i+1;
            c2 = j;
        elseif (theta >= -7*pi/8 && theta < -5*pi/8) || (theta >= 1*pi/8 && theta < 3*pi/8)
            r1 = i-1;
            c1 = j-1;
            r2 = i+1;
            c2 = j+1;
        end
        
        if (mag(i,j) > mag(r1,c1) && mag(i,j) > mag(r2,c2))
            suppressed(i,j) = mag(i,j);
        end
    end
end

low_thresh = 20;
high_thresh = 40;

edges = zeros(size(suppressed));
edges(suppressed >= high_thresh) = 1;
while true
    [r, c] = find(edges & suppressed >= low_thresh & ~edges);
    if isempty(r) || isempty(c)
        break;
    end
    for k = 1:numel(r)
        edges(r(k), c(k)) = 1;
    end
end

figure;
subplot(2,2,1); imshow(im); title('Original Image');
subplot(2,2,2); imshow(edges); title('Hysteresis thresheld binary edge image');

