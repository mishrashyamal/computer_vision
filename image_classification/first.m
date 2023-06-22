clear all;
close all;

image_directory = './CarData/TrainImages';
files = dir(image_directory);

X = [];
Y = [];
img_names = {};
for file = files'
    if ~file.isdir
        image = imread([image_directory, '/', file.name]);
        image = im2gray(image);
        numBins = 256;
        feature_vector = imhist(image, numBins)';
        if strcmp(file.name(1:3), 'pos')
            Y(end+1,1) = 1;
            img_names{end+1,1} = file.name;
        else
            Y(end+1,1) = 0;
            img_names{end+1,1} = file.name;
        end
        X(end+1,:) = feature_vector;
    end
end

rng(0);
indices = randperm(size(X,1));
num_samples = floor(size(X,1)/3);
X = X(indices,:);
Y = Y(indices,:);
Xtrain = X(1:2*num_samples,:);
Ytrain = Y(1:2*num_samples,:);
Xvalid = X(2*num_samples+1:end,:);
Yvalid = Y(2*num_samples+1:end,:);
img_names = img_names(indices,:);

k = 5;
num_validation_samples = size(Xvalid, 1);
predictions = zeros(num_validation_samples, 1);

for i = 1:num_validation_samples
    similarities = zeros(size(Xtrain, 1), 1);
    for j = 1:size(Xtrain, 1)
        hist_intersection = sum(min(Xvalid(i, :), Xtrain(j, :))) / numBins;
        similarities(j) = hist_intersection;
    end    
    [~, sorted_indices] = sort(similarities, 'descend');
    k_indices = sorted_indices(1:k);
    num_car = sum(Ytrain(k_indices));
    num_not_car = k - num_car;
    if num_car > num_not_car
        predictions(i) = 1; 
    else
        predictions(i) = 0;  
    end
end

accuracy = sum(predictions == Yvalid) / numel(Yvalid) * 100;
fprintf('Accuracy: %.2f%%\n', accuracy);

valid_files = img_names(2*num_samples+1:end,:);
Yvalid = Y(2*num_samples+1:end,:);

car_true_index = find(predictions == 1 & Yvalid == 1, 1, 'first');
not_car_true_index = find(predictions == 0 & Yvalid == 0, 1, 'first');
car_false_index = find(predictions == 1 & Yvalid == 0, 1, 'first');
not_car_false_index = find(predictions == 0 & Yvalid == 1, 1, 'first');

case_indices = [car_true_index, not_car_true_index, car_false_index, not_car_false_index];
case_labels = {'Correctly Labeled as Car', 'Correctly Labeled as Not a Car', 'Incorrectly Labeled as Car', 'Incorrectly Labeled as Not a Car'};

figure;
for i = 1:numel(case_indices)
    subplot(2, 2, i);
    imshow(imread(fullfile(image_directory, valid_files{case_indices(i)})));
    title(case_labels{i});
end


function hist = imhist(im, numBins)
    hist = zeros(1, numBins);
    [height, width, ~] = size(im);
    numPixels = height * width;

    for k = 1:numBins
        hist(k) = sum(sum(im == (k - 1))) / numPixels;
    end
end

