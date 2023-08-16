%% Part 2
inputImage = imread('cat.jpg');

h = size(inputImage, 1);
w = size(inputImage, 2);

% Change These Values for height and width
h_new = 330;
w_new = 250;

outputNN = uint8(zeros(h_new, w_new, 3));
outputBI = uint8(zeros(h_new, w_new, 3));

for y_new = 1:h_new
    for x_new = 1:w_new
        x = round(x_new * (w / w_new));
        y = round(y_new * (h / h_new));
        outputNN(y_new, x_new, :) = inputImage(y, x, :);
    end
end

for y_new = 1:h_new
    for x_new = 1:w_new
        x = x_new * (w / w_new);
        y = y_new * (h / h_new);
        
        x1 = floor(x);
        y1 = floor(y);
        x2 = ceil(x);
        y2 = ceil(y);
        
        dx = x - x1;
        dy = y - y1;
        
        top_interpolated = (1 - dx) * double(inputImage(y1, x1, :)) + dx * double(inputImage(y1, x2, :));
        bottom_interpolated = (1 - dx) * double(inputImage(y2, x1, :)) + dx * double(inputImage(y2, x2, :));
        
        interpolatedPixel = (1 - dy) * top_interpolated + dy * bottom_interpolated;
        
        outputBI(y_new, x_new, :) = uint8(interpolatedPixel);
    end
end

figure;
imshow(outputNN); title("Nearest Neighbor");

figure;
imshow(outputBI); title("BiLinear Interpolation");

%% Part 3

grayImage = rgb2gray(inputImage);

sigma = 1.5;

smoothedGrayImage = imgaussfilt(grayImage, sigma);

[Gx, Gy] = imgradientxy(smoothedGrayImage);

energy = sqrt(Gx.^2 + Gy.^2);
figure;
imshow(energy, []); title('Energy Function');

%% Part 4

% cost matrix
[rows, cols] = size(energy);
seamMatrix = zeros(rows, cols);
seamMatrix(1, :) = energy(1, :);
for i = 2:rows
    for j = 1:cols
        if j == 1
            seamMatrix(i, j) = energy(i, j) + min(seamMatrix(i-1, j), seamMatrix(i-1, j+1));
        elseif j == cols
            seamMatrix(i, j) = energy(i, j) + min(seamMatrix(i-1, j-1), seamMatrix(i-1, j));
        else
            seamMatrix(i, j) = energy(i, j) + min([seamMatrix(i-1, j-1), seamMatrix(i-1, j), seamMatrix(i-1, j+1)]);
        end
    end
end

%  optimal seam
[~, startCol] = min(seamMatrix(rows, :));
seamPath = zeros(rows, 1);
seamPath(rows) = startCol;
for i = rows-1:-1:1
    if seamPath(i+1) == 1
        [~, minIndex] = min([seamMatrix(i, 1), seamMatrix(i, 2)]);
        seamPath(i) = minIndex;
    elseif seamPath(i+1) == cols
        [~, minIndex] = min([seamMatrix(i, cols-1), seamMatrix(i, cols)]);
        seamPath(i) = cols - 1;
    else
        [~, minIndex] = min([seamMatrix(i, seamPath(i+1)-1), seamMatrix(i, seamPath(i+1)), seamMatrix(i, seamPath(i+1)+1)]);
        seamPath(i) = seamPath(i+1) - 2 + minIndex;
    end
end

% red line
seamThickness = 5;
overlayedImage = inputImage;
for i = 1:rows
    startCol = seamPath(i);
    for offset = -floor(seamThickness/2):floor(seamThickness/2)
        col = startCol + offset;
        col = max(1, min(col, cols)); 
        overlayedImage(i, col, 1) = 255;
        overlayedImage(i, col, 2) = 0;
        overlayedImage(i, col, 3) = 0;
    end
end

figure, imshow(overlayedImage);


%% Part 5

grayImage = rgb2gray(inputImage);

sigma = 1.5;

smoothedGrayImage = imgaussfilt(grayImage, sigma);

[Gx, Gy] = imgradientxy(smoothedGrayImage);

energy = sqrt(Gx.^2 + Gy.^2);

% cost matrix
[rows, cols] = size(energy);
seamMatrix = zeros(rows, cols);
seamMatrix(1, :) = energy(1, :);
for i = 2:rows
    for j = 1:cols
        if j == 1
            seamMatrix(i, j) = energy(i, j) + min(seamMatrix(i-1, j), seamMatrix(i-1, j+1));
        elseif j == cols
            seamMatrix(i, j) = energy(i, j) + min(seamMatrix(i-1, j-1), seamMatrix(i-1, j));
        else
            seamMatrix(i, j) = energy(i, j) + min([seamMatrix(i-1, j-1), seamMatrix(i-1, j), seamMatrix(i-1, j+1)]);
        end
    end
end
% Change Below
outputVideo = VideoWriter('bikeSeam.mp4', 'MPEG-4');
open(outputVideo);
targetSize = [size(inputImage, 1), size(inputImage, 2)];
while cols > 1
    %  optimal seam
    [~, startCol] = min(seamMatrix(rows, :));
    seamPath = zeros(rows, 1);
    seamPath(rows) = startCol;
    for i = rows-1:-1:1
        prevCol = seamPath(i+1); 
        if prevCol == 1
            [~, minIndex] = min([seamMatrix(i, 1), seamMatrix(i, 2)]);
            seamPath(i) = minIndex;
        elseif prevCol == cols
            [~, minIndex] = min([seamMatrix(i, cols-1), seamMatrix(i, cols)]);
            seamPath(i) = cols - 1;
        else
            [~, minIndex] = min([seamMatrix(i, prevCol - 1), seamMatrix(i, prevCol), seamMatrix(i, prevCol + 1)]);
            seamPath(i) = prevCol - 2 + minIndex;
        end
    end
    
    seamThickness = 5; 
    overlayedImage = inputImage;
    for i = 1:rows
        startCol = seamPath(i);
        for offset = -floor(seamThickness/2):floor(seamThickness/2)
            col = startCol + offset;
            col = max(1, min(col, cols)); 
            overlayedImage(i, col, 1) = 255;
            overlayedImage(i, col, 2) = 0;
            overlayedImage(i, col, 3) = 0;
        end
    end
    
    padAmount = targetSize(2) - size(overlayedImage, 2);
    paddedOverlayedImage = padarray(overlayedImage, [0, padAmount], 0, 'post');
    
    writeVideo(outputVideo, paddedOverlayedImage);
    
    inputImage(:, startCol, :) = [];
    energy(:, startCol) = [];
    seamMatrix(:, startCol) = [];
    
    cols = cols - 1;
end

close(outputVideo);

