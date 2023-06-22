clear all;

image = imread('test.png');

if size(image, 3) == 3
    image = rgb2gray(image);
end

kernel_size = 8;
sigma = 1.5;

low_threshold = 0.2 ;
high_threshold = 0.4;

% find edges of the image
edges = getEdgesImage(image, kernel_size, sigma, low_threshold, high_threshold);

figure;
imshow(edges)
title('edges of the image');

% calcuate Hough Transform for Line Detection 
theta_steps = 180;
hough_hist = calculateHoughTransform(edges, theta_steps);

figure;
imshow(hough_hist, [])
title('Hough Transform for Line Detection');

% filter relevent lines
[numLines, lines] = indentifyReleventLine(hough_hist);
maxima_cols = lines(:, 1);
maxima_rows = lines(:, 2);

% plot the maxima
figure;
imshow(hough_hist, [])
title('Hough Transform with Local Maxima');

hold on;
plot(maxima_cols, maxima_rows, 'r.', 'MarkerSize', 5);
hold off;

% superimpose the potential lines 
thresholdDistance = 20;
selectedLines = superimposeSelectedLines(edges, lines, thresholdDistance);
rho_steps = round(sqrt(size(edges, 1)^2 + size(edges, 2)^2));

figure;
imshow(edges);

for k = 1:size(selectedLines, 1)
    theta_rad = selectedLines(k, 2) * pi / 180;
    rho = selectedLines(k, 1) - rho_steps - 1;

    x = 1:size(edges, 2);
    y = (rho - x*cos(theta_rad)) / sin(theta_rad);

    line(x, y, 'Color', 'r', 'LineWidth', 1);
end


% find intersection points
intersectionPoints = findIntersectionPoints(selectedLines,edges);

% group intersection points which is required for image rectification
finalintersectionPoints = calculateFinalIntersectionPoints(intersectionPoints, image);

% image rectfication
outputWidth = 8.5 * 100;  
outputHeight = 11 * 100; 
rectifiedImage = rectifyImage(image, finalintersectionPoints, outputWidth, outputHeight);

figure;
imshow(flipud(rectifiedImage));

function rectifiedImage = rectifyImage(image, finalintersectionPoints, outputWidth, outputHeight)
    rectifiedImage = zeros(outputHeight, outputWidth, size(image, 3), 'uint8');

    p1 = finalintersectionPoints(1,:);
    p2 = finalintersectionPoints(2,:);
    p3 = finalintersectionPoints(3,:);
    p4 = finalintersectionPoints(4,:);

    q1 = [1, 1];
    q2 = [outputWidth, 1];
    q3 = [outputWidth, outputHeight];
    q4 = [1, outputHeight];

    H = computeHomography([p1; p2; p3; p4], [q1; q2; q3; q4]);

    for y = 1:outputHeight
        for x = 1:outputWidth
            p = H \ [x; y; 1];
            p = p / p(3);

            if p(1) >= 1 && p(1) <= size(image, 2) && p(2) >= 1 && p(2) <= size(image, 1)
                rectifiedImage(y, x, :) = image(round(p(2)), round(p(1)), :);
            end
        end
    end
end

function H = computeHomography(p, q)
    numPoints = size(p, 1);

    A = zeros(numPoints * 2, 9);
    for i = 1:numPoints
        x = p(i, 1);
        y = p(i, 2);
        u = q(i, 1);
        v = q(i, 2);

        A(2*i-1, :) = [-x, -y, -1, 0, 0, 0, u*x, u*y, u];
        A(2*i, :) = [0, 0, 0, -x, -y, -1, v*x, v*y, v];
    end

    [~, ~, V] = svd(A);
    h = V(:, end);
    H = reshape(h, 3, 3)';
end


function [slope, intercept] = extractSlopeIntercept(equationVariable)
    equation = char(equationVariable);
    pattern = '[-+]?\d*\.?\d+';
    matches = regexp(equation, pattern, 'match');
    slope = NaN;
    intercept = NaN;
   
    if numel(matches) == 2
        slope = str2double(matches{1});
        intercept = str2double(matches{2});
       
        if isnan(slope) || isnan(intercept)
            slope = NaN;
            intercept = NaN;
        end
    end
end


function [numLines, lines] = indentifyReleventLine(houghTransformLine)
    threshold = 0.5 * max(houghTransformLine(:));
    [h, w] = size(houghTransformLine);
    localMaxima = findLocalMaxima(houghTransformLine);
    linesAboveThreshold = houghTransformLine >= threshold;
    [rows, cols] = find(localMaxima & (linesAboveThreshold));
    [~, indices] = sort(houghTransformLine(sub2ind([h, w], rows, cols)), 'descend');
    rows = rows(indices);
    cols = cols(indices);
    numLines =  numel(indices);
    lines = [cols(1:numLines), rows(1:numLines)];
end

function localMaxima = findLocalMaxima(houghTransform)
    [rows, cols] = size(houghTransform);
    localMaxima = false(size(houghTransform));

    for i = 2:rows-1
        for j = 2:cols-1
            pixelValue = houghTransform(i, j);
            if pixelValue > houghTransform(i-1, j-1) && ...
               pixelValue > houghTransform(i-1, j) && ...
               pixelValue > houghTransform(i-1, j+1) && ...
               pixelValue > houghTransform(i, j-1) && ...
               pixelValue > houghTransform(i, j+1) && ...
               pixelValue > houghTransform(i+1, j-1) && ...
               pixelValue > houghTransform(i+1, j) && ...
               pixelValue > houghTransform(i+1, j+1)
                localMaxima(i, j) = true;
            end
        end
    end
end


function [labels, numLabels] = bwlabel(binary_image)
    [height, width] = size(binary_image);
    labels = zeros(height, width);
    current_label = 1;
    
    for i = 1:height
        for j = 1:width
            if binary_image(i, j) == 1 && labels(i, j) == 0
                stack = [i, j];
                while ~isempty(stack)
                    current_pixel = stack(1, :);
                    stack(1, :) = [];
                    row = current_pixel(1);
                    col = current_pixel(2);

                    if row >= 1 && row <= height && col >= 1 && col <= width && binary_image(row, col) == 1 && labels(row, col) == 0
                        labels(row, col) = current_label;
                        stack = [stack; row-1, col; row+1, col; row, col-1; row, col+1];
                    end
                end
                current_label = current_label + 1;
            end
        end
    end
   
    numLabels = current_label - 1;
end

function areaStats = calculateRegionArea(labels)
    [height, width] = size(labels);
   
    uniqueLabels = unique(labels(:));
    numLabels = numel(uniqueLabels) - 1;  
    areaStats = struct('Area', zeros(numLabels, 1), 'Label', zeros(numLabels, 1));
   
    for i = 1:numLabels
        currentLabel = uniqueLabels(i+1);  
        area = sum(labels(:) == currentLabel);
        areaStats(i).Area = area;
        areaStats(i).Label = currentLabel;
    end
end

function edges = getEdgesImage(image, kernel_size, sigma, low_threshold, high_threshold)
    if size(image, 3) == 3
        image = rgb2gray(image);
    end

    radius = floor(kernel_size / 2);
    [X, Y] = meshgrid(-radius:radius, -radius:radius);
    gaussian_kernel = exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
    smoothed_image = conv2(double(image), double(gaussian_kernel), 'same');

    kernel_x = [-1 0 1; -2 0 2; -1 0 1];
    kernel_y = [1 2 1; 0 0 0; -1 -2 -1];
    Gx = conv2(smoothed_image, kernel_x, 'same');
    Gy = conv2(smoothed_image, kernel_y, 'same');

    g_mag = sqrt(Gx.^2 + Gy.^2);
    g_orien = zeros(size(Gx));
   
    g_orien(Gx == 0 & Gy > 0) = pi / 2;
    g_orien(Gx == 0 & Gy < 0) = -pi / 2;
    g_orien(Gy == 0 & Gx > 0) = 0;
    g_orien(Gy == 0 & Gx < 0) = pi;
    g_orien(Gx > 0 & Gy > 0) = atan(Gy(Gx > 0 & Gy > 0) ./ Gx(Gx > 0 & Gy > 0));
    g_orien(Gx < 0 & Gy > 0) = pi + atan(Gy(Gx < 0 & Gy > 0) ./ Gx(Gx < 0 & Gy > 0));
    g_orien(Gx < 0 & Gy < 0) = -pi + atan(Gy(Gx < 0 & Gy < 0) ./ Gx(Gx < 0 & Gy < 0));
    g_orien(Gx > 0 & Gy < 0) = atan(Gy(Gx > 0 & Gy < 0) ./ Gx(Gx > 0 & Gy < 0));

    suppressed = zeros(size(g_mag));
    for i = 2:size(g_mag, 1)-1
        for j = 2:size(g_mag, 2)-1
            theta = g_orien(i,j);
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

            if (g_mag(i,j) > g_mag(r1,c1) && g_mag(i,j) > g_mag(r2,c2))
                suppressed(i,j) = g_mag(i,j);
            end
        end
    end

    low_threshold = low_threshold * max(suppressed(:));
    high_threshold = high_threshold * max(suppressed(:));
    edges = zeros(size(suppressed));

    for i = 2:size(suppressed, 1)-1
        for j = 2:size(suppressed, 2)-1
            if g_mag(i,j) > high_threshold
                edges(i,j) = 1;
            elseif g_mag(i,j) > low_threshold && ...
                    (g_mag(i-1,j-1) > high_threshold || ...
                    g_mag(i-1,j) > high_threshold || ...
                    g_mag(i-1,j+1) > high_threshold || ...
                    g_mag(i,j-1) > high_threshold || ...
                    g_mag(i,j+1) > high_threshold || ...
                    g_mag(i+1,j-1) > high_threshold || ...
                    g_mag(i+1,j) > high_threshold || ...
                    g_mag(i+1,j+1) > high_threshold)
                edges(i,j) = 1;
            end
        end
    end
end




function hough_hist = calculateHoughTransform(edges, theta_steps)

    [labels, numLabels] = bwlabel(edges);
    edgeStats = calculateRegionArea(labels);
    [~, idx] = max([edgeStats.Area]);
    edges = (labels == idx);
    rho_steps = round(sqrt(size(edges, 1)^2 + size(edges, 2)^2));
    hough_hist = zeros(theta_steps, 2*rho_steps + 1);
    
    for y = 1:size(edges, 1)
        for x = 1:size(edges, 2)
            if edges(y, x) == 1
                for theta_deg = 1:theta_steps
                    theta_rad = (theta_deg - 1) * pi / theta_steps;
                    rho = round(x * cos(theta_rad) + y * sin(theta_rad));
                    hough_hist(theta_deg, rho + rho_steps + 1) = hough_hist(theta_deg, rho + rho_steps + 1) + 1;
                end
            end
        end
    end
end


function selectedLines = superimposeSelectedLines(edges, lines, thresholdDistance)
    selectedLines = [];
    numLines = size(lines, 1);

    for i = 1:numLines
        line = lines(i, :);
        distances = [];
        for j = 1:size(selectedLines, 1)
            selectedLine = selectedLines(j, :);
            distance = calculateDistance(line, selectedLine);
            distances = [distances, distance];
        end

        if isempty(selectedLines) || all(distances > thresholdDistance)
            selectedLines = [selectedLines; line];
        end
    end
end

function distance = calculateDistance(line1, line2)
    distance = norm(line1 - line2);
end


function intersectionPoints = findIntersectionPoints(selectedLines,edges)
    imageSize = size(image);
    equations = cell(size(selectedLines, 1), 1);
    rho_steps = round(sqrt(size(edges, 1)^2 + size(edges, 2)^2));
 
    for k = 1:size(selectedLines, 1)
        theta_rad = selectedLines(k, 2) * pi / 180;
        rho = selectedLines(k, 1) - rho_steps - 1;
        m = -cos(theta_rad) / sin(theta_rad);
        b = rho / sin(theta_rad);
        equations{k} = sprintf('y = %.2f * x + %.2f', m, b);
        x = 1:imageSize(2);
        y = m * x + b;
    end

    intersectionPoints = [];

    for i = 1:numel(equations)-1
        equation1 = equations{i};
        [m1, b1] = extractSlopeIntercept(equation1); 

        for j = i+1:numel(equations)
            equation2 = equations{j};
            [m2, b2] = extractSlopeIntercept(equation2); 

            x_intersect = (b2 - b1) / (m1 - m2);
            y_intersect = m1 * x_intersect + b1;

            intersectionPoints = [intersectionPoints; x_intersect, y_intersect];
        end
    end
end


function finalintersectionPoints = calculateFinalIntersectionPoints(intersectionPoints, image)
    sortedIntersectionPoints = sortrows(intersectionPoints, 1);
    imageWidth = size(image, 2);
    imageHeight = size(image, 1);

    leftmostIntersectionPoints = sortedIntersectionPoints(1:3, :);
    rightmostIntersectionPoints = sortedIntersectionPoints(end-2:end, :);

    finalintersectionPoints = [];

    for i = 1:size(leftmostIntersectionPoints, 1)
        x1 = leftmostIntersectionPoints(i, 1);
        x = leftmostIntersectionPoints(i, 1) / imageWidth;
        y1 = leftmostIntersectionPoints(i, 2);
        y = leftmostIntersectionPoints(i, 2) / imageHeight;

        if x >= 0 && x <= imageWidth && y >= 0 && y <= imageHeight
            finalintersectionPoints = [finalintersectionPoints; x1, y1];
        end
    end

    for i = 1:size(rightmostIntersectionPoints, 1)
        x1 = rightmostIntersectionPoints(i, 1);
        x = rightmostIntersectionPoints(i, 1) / imageWidth;
        y1 = rightmostIntersectionPoints(i, 2);
        y = rightmostIntersectionPoints(i, 2) / imageHeight;

        if x >= 0 && x <= imageWidth && y >= 0 && y <= imageHeight
            finalintersectionPoints = [finalintersectionPoints; x1, y1];
        end
    end

    center = mean(finalintersectionPoints);
    angles = atan2d(finalintersectionPoints(:,2)-center(2), finalintersectionPoints(:,1)-center(1));
    [~, order] = sort(angles, 'descend');
    sortedPoints = finalintersectionPoints(order, :);
    finalintersectionPoints = sortedPoints;
end


