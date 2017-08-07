%Load in images

imnames = {'atrium/IMG_1347.JPG','atrium/IMG_1348.JPG','atrium/IMG_1349.JPG'};
nimages = length(imnames);
baseim = 1; %index of the central "base" image

for i = 1:nimages
  ims{i} = imresize(im2double(imread(imnames{i})), 0.25);
  ims_gray{i} = rgb2gray(ims{i});
  [h(i),w(i),~] = size(ims{i});
end

% get corresponding points between each image and the central base image
% 
% for i = 1:nimages
%    if (i ~= baseim)
%      %run interactive select tool to click corresponding points on base and non-base image
%      [moving_points(:, :, i - 1), fixed_points(:, :, i - 1)] = cpselect(ims{i}, ims{baseim}, 'Wait', true);
% 
%      %refine the user clicks using cpcorr
%       moving_points(:, :, i - 1) = cpcorr(moving_points(:, :, i - 1), fixed_points(:, :, i - 1), ims{i}(:, :, 1), ims{i}(:, :, 1));
%    end
% end

%save atriumpts.mat moving_points fixed_points
%load atriumpts.mat

for i = 1:nimages
    if (i ~= baseim)
        H(:, :, i) = computeHomography(fixed_points(:, 1, i - 1)', fixed_points(:, 2, i - 1)', moving_points(:, 1, i - 1)', moving_points(:, 2, i - 1)');
    else
        H(:, :, i) = eye(3)
    end
end


% compute where corners of each warped image end up
xmin = Inf;
ymin = Inf;
xmax = -Inf;
ymax = -Inf;

for i = 1:nimages
    cx = [1, 1, w(i), w(i)];
    cy = [1, h(i), 1, h(i)];
    [cx_warped, cy_warped] = applyHomography(H(:, :, i), cx, cy);
     
% find corners of a rectangle that contains all the warped image
%  corner points
    xmin = floor(min([cx_warped, xmin]));
    xmax = ceil(max([cx_warped, xmax]));
    ymin = floor(min([cy_warped, ymin]));
    ymax = ceil(max([cy_warped, ymax])); 
end

% Use H and interp2 to perform inverse-warping of the source image to align it with the base image

[xx, yy] = meshgrid(xmin:xmax, ymin:ymax);
[pp, qq] = meshgrid(1:xmax-xmin, 1:ymax-ymin);

for i = 1:nimages
    [xq,yq] = applyHomography(inv(H(:,:,i)), xx(:)', yy(:)');
    R = interp2(ims{i}(:, :, 1), xq, yq);
    R = reshape(R, size(xx));
    G = interp2(ims{i}(:, :, 2), xq, yq);
    G = reshape(G, size(xx));
    B = interp2(ims{i}(:, :, 3), xq, yq);
    B = reshape(B, size(xx));
    J{i} = cat(3, R, G, B);
    mask{i} = ~isnan(R);  %interp2 puts NaNs outside the support of the warped image
    J{i}(isnan(J{i})) = 0;
end

mask1and2 = mask{1} & mask{2};
mask1and3 = mask{1} & mask{3};
mask2and3 = mask{2} & mask{3};

beforeblend1= J{1};
beforeblend2 = J{2};
beforeblend3 = J{3};

%put these in writeup along with blended image
figure,imshow(beforeblend1)
figure,imshow(beforeblend2)
figure,imshow(beforeblend3)

%Create the gaussian filter with hsize = [450 450] and sigma = 160
gausFilter = fspecial('gaussian', [450 450], 160);

for i = 1:nimages
   % blur and clip mask{i} to get an alpha map for each image
   alpha{i} = imfilter(im2double(mask{i}), gausFilter);
   alpha{i}(~mask{i}) = 0;
end

% only blur where images intersect
alpha{1}(alpha{1} & ~mask1and2 & ~mask1and3) = 1;
alpha{2}(alpha{2} & ~mask1and2 & ~mask2and3) = 1;
alpha{3}(alpha{3} & ~mask1and3 & ~mask2and3) = 1;
 

A = alpha{1} + alpha{2} + alpha{3};
pixels_to_fix = find(A ~= 0 & A ~= 1);

% scale alpha maps to sum to 1 at every pixel location
for i = 1:size(pixels_to_fix)
    index = pixels_to_fix(i);
    sum = alpha{1}(index) + alpha{2}(index) + alpha{3}(index);
    alpha{1}(index) = (alpha{1}(index) / sum);
    alpha{2}(index) = (alpha{2}(index) / sum);
    alpha{3}(index) = (1 - (alpha{1}(index) + alpha{2}(index)) );
end

% finally blend together the resulting images into the final mosaic
K = zeros(size(J{1}));
for i = 1:nimages
    K = K + repmat(alpha{i}, [1, 1, 3]).*J{i};
end

% display the result
figure(1), 
imagesc(K); axis image;
imshow(K)