clear all; close all; clc
warning('off', 'images:initSize:adjustingMag')
%% Section 1Formating and Inputs
close all

calib_img = imread(['.' filesep 'calibration_image' filesep 'DSC_0001.JPG']); % Read image
h1=figure(1); imshow(calib_img,'InitialMagnification',40); % Show Entire Image
[y_max_pxl, x_max_pxl] = size(calib_img(:,:,2)); % find size of image

% INPUT VALUES %
rotation_angle = 0.6; % rotation of image
number_of_calibration_points = 20; % number of points on ruler
limval=20; %crop horizontally of image
crop_box = [limval 875 x_max_pxl-limval 250]; %% [XMIN YMIN WIDTH HEIGHT];

% Rotate and Calibrate
calib_img = imrotate(calib_img, rotation_angle, 'loose', 'bilinear'); % Rotates image
calib_img = imadjust(calib_img, stretchlim(calib_img),[]); % adusts contras

% Crop
calib_img = calib_img(:,:,2); % Choosen the second channel of the RGB Image
calib_img = imcrop(calib_img,crop_box); % Crop Image

imshow(calib_img,'InitialMagnification',120); hold on % Show rotate, cropped and adjusted image
plot([0 x_max_pxl], [200 200]) % plot line for reference
hold off

%% Section 2 Calibration points on image and in real world
%close all

figure; imshow(calib_img,'InitialMagnification',60); hold on % show image

world_x = linspace(0,19,number_of_calibration_points)'*10; % Create world scale 0-190 cm

[image_x,~] = ginput(number_of_calibration_points); % graphically input corrdinates

format long % show extra decimals for when saving the coordinates
%% Section 3 display chosen calibration points

format compact % make displaying format compact again to make it manageale

figure; imshow(calib_img,'InitialMagnification',60); hold on % show figure

world_x = linspace(0,19,number_of_calibration_points)'*10; % Show real world scale

%This is where you save your pixel coordinates from Section 2
image_x =[0.239468202668891
   0.369645433694745
   0.507261363636364
   0.642397727272727
   0.781253440366973
   0.921348936613845
   1.065163782318599
   1.206499061718098
   1.351553690575479
   1.495368536280234
   1.640423165137615
   1.785477793994996
   1.928052856547123
   2.070627919099250
   2.209483632193495
   2.348339345287740
   2.484715492076731
   2.619851855713095
   2.752508653044204
   2.886405233527941
]*1e3;

for idx = 1: numel(image_x) % plotting the points on the image
  plot([image_x(idx) image_x(idx)], [0 size(calib_img,2)])
end
hold off



%% calibration fit (adjust for distortion in x)
fit_x = fit(image_x, world_x, 'poly2'); % make fit of image pixels to real world coordinates

figure; % Make figure
plot(fit_x, image_x, world_x)
xlabel('x[pxl]');
ylabel('x[cm]')