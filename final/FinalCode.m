%% Part 1
clear; clc;
image = imread('final\paper_reduced.jpg');
grayImage = rgb2gray(image);

kernel = (1/9) * ones(3, 3); 
edgeX = [-1, 0, 1; -2, 0, 2; -1, 0, 1];  
edgeY = [-1, -2, -1; 0, 0, 0; 1, 2, 1];  

smoothed = conv2(double(grayImage), kernel, 'same');
gradX = conv2(smoothed, edgeX, 'same');
gradY = conv2(smoothed, edgeY, 'same');
gradM = sqrt(gradX.^2 + gradY.^2);

threshold = 200; %200 for blank page, 300 for notes page
edgeMap = gradM > threshold;
figure;
imshow(edgeMap); title('Binary Edge Map');
%% Part 2
theta_step = 1;
rho_step = 1;
maxRho = round(hypot(size(grayImage, 1), size(grayImage, 2)));
thetaRange = 0:theta_step:180;
rhoRange = -maxRho:rho_step:maxRho;
hough = zeros(length(thetaRange), length(rhoRange));
[y_coords, x_coords] = find(edgeMap);
num_points = length(x_coords);

for i = 1:num_points
    x = x_coords(i);
    y = y_coords(i);
    
    for thetaIndex = 1:length(thetaRange)
        theta = deg2rad(thetaRange(thetaIndex));
        rho = round(x * cos(theta) + y * sin(theta));
        rhoIndex = round((rho + maxRho) / rho_step) + 1;
        hough(thetaIndex, rhoIndex) = hough(thetaIndex, rhoIndex) + 1;
    end
end

figure;
imagesc(rhoRange, thetaRange, hough);
colormap(gca, gray);
%% Part 3

potential = zeros(0, 2); 
for thetaIndex = 1:length(thetaRange)
    for rhoIndex = 1:length(rhoRange)
        if hough(thetaIndex, rhoIndex) > threshold
            isMax = true;
            for i = -1:1
                for j = -1:1
                    if (i ~= 0 || j ~= 0) && ...
                       (rhoIndex + i >= 1 && rhoIndex + i <= length(rhoRange)) && ...
                       (thetaIndex + j >= 1 && thetaIndex + j <= length(thetaRange))
                        if hough(thetaIndex, rhoIndex) <= hough(thetaIndex + j, rhoIndex + i)
                            isMax = false;
                            break;
                        end
                    end
                end
                if ~isMax
                    break;
                end
            end
            
            if isMax
                potential = [potential; thetaRange(thetaIndex), rhoRange(rhoIndex)];
            end
        end
    end
end


figure;
imshow(edgeMap);
hold on;

for i = 1:size(potential, 1)
    theta = deg2rad(potential(i, 1));
    rho = potential(i, 2);
    x1 = 0;
    x2 = size(image, 2);
    y1 = (rho - x1 * cos(theta)) / sin(theta);
    y2 = (rho - x2 * cos(theta)) / sin(theta);

    plot([x1, x2], [y1, y2], 'r', 'LineWidth', 1);
end

hold off;

%% part 4
intersections = [];

for i = 1:size(potential, 1)
    for j = i+1:size(potential, 1)
        theta1 = deg2rad(potential(i, 1));
        rho1 = potential(i, 2);
        theta2 = deg2rad(potential(j, 1));
        rho2 = potential(j, 2);
        A = [cos(theta1), sin(theta1); cos(theta2), sin(theta2)];
        b = [rho1; rho2];
        intersect = A \ b;
        xIntersect = intersect(1);
        yIntersect = intersect(2);
        if xIntersect >= 1 && xIntersect <= size(image, 2) && yIntersect >= 1 && yIntersect <= size(image, 1)
            intersections = [intersections; xIntersect, yIntersect];
        end
    end
end

if ~isempty(intersections)
    sorted = sortrows(intersections, 1);
    corner1 = sorted(1, :);
    corner3 = sorted(end, :);
    sorted = sortrows(intersections, 2);
    corner2 = sorted(1, :);
    corner4 = sorted(end, :);
    figure;
    imshow(image)
    hold on;
    plot(corner1(1), corner1(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(corner3(1), corner3(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(corner2(1), corner2(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot(corner4(1), corner4(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    hold off;
    
end
%% part 5
originalPoints = [corner1;corner2;corner3; corner4];
[heightPixels, widthPixels, numberOfColorChannels] = size(image);

blankImage = uint8(zeros(heightPixels, widthPixels, 3));

blankPoints = [0, 0; widthPixels, 0; widthPixels, heightPixels; 0, heightPixels];

H = calculateHomographyMatrix(blankPoints, originalPoints);

for y = 1:heightPixels
    for x = 1:widthPixels
        originalLocation = H * [x; y; 1];
        originalX = originalLocation(1) / originalLocation(3);
        originalY = originalLocation(2) / originalLocation(3);
        
        if originalX >= 1 && originalX <= size(image, 2) && originalY >= 1 && originalY <= size(image, 1)
            blankImage(y, x, :) = image(round(originalY), round(originalX), :);
        end
    end
end

figure
imshow(blankImage);




function H = calculateHomographyMatrix(points1, points2)
    numPoints = size(points1, 1);
    A = zeros(2 * numPoints, 9);

    for i = 1:numPoints
        X1 = points1(i, 1);
        Y1 = points1(i, 2);
        X2 = points2(i, 1);
        Y2 = points2(i, 2);

        A(2 * i - 1, :) = [-X1, -Y1, -1, 0, 0, 0, X1 * X2, Y1 * X2, X2];
        A(2 * i, :) = [0, 0, 0, -X1, -Y1, -1, X1 * Y2, Y1 * Y2, Y2];
    end
    [~, ~, V] = svd(A);
    h = V(:, 9);
    H = reshape(h, [3, 3])';
end

