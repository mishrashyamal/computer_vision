
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

imshow(my_img);
title('Binary Image with Line and Circle');




