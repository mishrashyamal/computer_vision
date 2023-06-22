
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


[height, width] = size(my_img);
rMax = round(sqrt(height^2 + width^2));
hough_line = zeros(180, rMax*2);
max = 0;

for i = 1:width
    for j = 1:height
        if my_img(j,i) == 1
            for deg = 1:180
                theta = (deg/180) * pi;
                r = round(cos(theta)*i + sin(theta)*j);
                hough_line(deg, r+rMax) = hough_line(deg, r+rMax) + 1;
                if hough_line(deg, r+rMax) > max
                    max = hough_line(deg, r+rMax);
                end
            end
        end
    end
end

for deg = 1:180
    for r = 1:rMax*2
        hough_line(deg, r) = hough_line(deg, r)/max;
    end
end


figure; 
imshow(hough_line, []);
title('Hough Transform for a Line');








