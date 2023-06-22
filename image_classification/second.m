close;
clear all;

image_directory = './CarData/TrainImages';
files = dir(image_directory);

X = [];
Y = [];
img_names = {};
for f=files'
    if ~f.isdir
        im = imread([image_directory, '/', f.name]);
        im = im2gray(im);

        fv = [];
        subImageSize = [20,20];
        numberSubImage = 10;
        for x = 1:numberSubImage:size(im,1) - subImageSize(1)+1
            for y = 1:numberSubImage:size(im,2) - subImageSize(2)+1
                subImage = im(x:x+subImageSize(1)-1, y:y+subImageSize(2)-1);
                finalHOG = computeHOG(subImage,8);

                fv = [fv,finalHOG];
            end
        end
        
        if strcmp(f.name(1:3),'pos')
            Y(end+1,1) = 1;
            img_names{end+1,1} = f.name;
        else 
            Y(end+1,1) = 0;
            img_names{end+1,1} = f.name;
        end
        X(end+1, :) = fv(:)';
    end
end


rng(0);
inds = randperm(size(X,1));
num = floor(size(X,1)/3);
X = X(inds,:);
Y = Y(inds,:);
Xtrain = X(1:2*num,:);
Ytrain = Y(1:2*num,:);
Xvalid = X(2*num+1:end,:);
Yvalid = Y(2*num+1:end,:);
img_names = img_names(inds,:);

k = 5;
prediction = zeros(size(Yvalid));

for i=1:size(Xvalid,1)
    a = Xvalid(i,:);
    distance = zeros(size(Xtrain,1),1);

    for j = 1:size(Xtrain,1)
        b = Xtrain(j,:);
        distance(j) = sum(min(a,b));
    end
     [~,index] = sort(distance,'descend');
     knn_labels = Ytrain(index(1:k));
     prediction(i) = mode(knn_labels);
end

accuracy = sum(prediction == Yvalid) / numel(Yvalid) * 100;
fprintf('Accuracy: %.2f%%\n', accuracy);

Valid_img = img_names(2*num+1:end,:);
Yvalid = Y(2*num+1:end,:);


car_true_index = find(prediction==1&Yvalid==1,1);
not_car_true_index = find(prediction==0&Yvalid==0,1);
car_false_index = find(prediction==1&Yvalid==0,1);
not_car_false_index = find(prediction==0&Yvalid==1 ,1);
 
figure;
subplot(2, 2, 1);
imshow(imread(fullfile(image_directory, Valid_img{car_true_index})));
title('Correctly Labeled as Car');
subplot(2, 2, 2);
imshow(imread(fullfile(image_directory, Valid_img{not_car_true_index})));
title('Correctly Labeled as Not a Car');
subplot(2, 2, 3);
imshow(imread(fullfile(image_directory, Valid_img{car_false_index})));
title('Incorrectly Labeled as Car');
subplot(2, 2, 4);
imshow(imread(fullfile(image_directory, Valid_img{not_car_false_index})));
title('Incorrectly Labeled as Not a Car');

 function hogFeatures = computeHOG(image, numBins)
    [Gmag, Gdir] = imgradient(image);
    
    cellSize = [8, 8];
    [numRows, numCols] = size(Gmag);
    numCellsRow = floor(numRows / cellSize(1));
    numCellsCol = floor(numCols / cellSize(2));
    
    hogFeatures = zeros(numCellsRow, numCellsCol, numBins);
    binEdges = linspace(0, 180, numBins + 1);
    
    for row = 1:numCellsRow
        for col = 1:numCellsCol
            cellMag = Gmag((row-1)*cellSize(1)+1:row*cellSize(1), ...
                           (col-1)*cellSize(2)+1:col*cellSize(2));
            cellDir = Gdir((row-1)*cellSize(1)+1:row*cellSize(1), ...
                           (col-1)*cellSize(2)+1:col*cellSize(2));
            
            histogram = histcounts(cellDir(:), binEdges, 'Normalization', 'probability');
            
            hogFeatures(row, col, :) = histogram;
        end
    end
    
    hogFeatures = reshape(hogFeatures, [], numBins);
    hogFeatures = hogFeatures(:)';
end