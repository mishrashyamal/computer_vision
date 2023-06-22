
my_img = zeros(400, 400);

for x = 1:400
    y = round(x - 100);
    if y > 0 && y <= 400
        my_img(y, x) = 1;
    end
end

for theta = 0:0.01:(2*pi)
    x = round(100 + 50*cos(theta));
    y = round(200 + 50*sin(theta));
    if x > 0 && x <= 400 && y > 0 && y <= 400
        my_img(y, x) = 1;
    end
end

rMin = 30;
rMax = 80;
numRadii = rMax - rMin + 1;
hough_circle = zeros(size(my_img, 1), size(my_img, 2), numRadii);

for r = rMin:rMax
    for y0 = 1:size(my_img,1)
        for x0 = 1:size(my_img,2)
            if my_img(y0,x0) == 1 
                for theta = 0:pi/180:2*pi 
                    x = round(x0 + r*cos(theta)); 
                    y = round(y0 + r*sin(theta)); 
                    if x > 0 && x <= size(my_img,2) && y > 0 && y <= size(my_img,1)
                        hough_circle(y,x,r-rMin+1) = hough_circle(y,x,r-rMin+1) + 1;
                    end
                end
            end
        end
    end
end


[maxVal, maxIdx] = max(hough_circle(:));
[x0Idx, y0Idx, rIdx] = ind2sub(size(hough_circle), maxIdx);
x0 = x0Idx;
y0 = y0Idx;
r = rIdx + rMin - 1;



figure;
imshow(hough_circle(:,:,rIdx), 'DisplayRange', [0.5 max(hough_circle(:))]);
title(sprintf('Hough Transform for r = %d', r));

fprintf('Values:\n');
fprintf('(x0, y0, r) = (%d, %d, %d)\n', x0, y0, r);



