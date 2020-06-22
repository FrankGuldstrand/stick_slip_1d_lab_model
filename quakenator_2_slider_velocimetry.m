%% Section 1
close all % close all figures
warning('off', 'images:initSize:adjustingMag')

number_of_sliders = 2; % How many sliding objects
image_folder = ['.' filesep 'experiment']; % path to experiment images

d = dir([image_folder, filesep, '*.JPG']); % checks the contents of JPGs in folder
number_of_time_steps = numel(d); % counts the JPGs in the folder

p= zeros(number_of_time_steps,number_of_sliders); % Variable for storing force gauge, 1 fg, 2b, 3b

%% Section 2 Adjusting and testing the settings of image analysis

% INPUTS
crop_box = [limval 910 x_max_pxl-limval 110]; % New Smaller crop box.  %% [XMIN YMIN WIDTH HEIGHT];
time_step = 25; % Time step to test ( from number_of_time_steps)
black_white_threshold = 0.3; %IMPORTANT BLACK WHITE THRESHOLD

  disp(['Processing: ' d(time_step).name ' ...']) % Display which files is being shown in command window
  
  img = imread([image_folder, filesep, d(time_step).name]); % Load Image  
  img = img(:,:,2); % Choose one channel of color
  
  img = imrotate(img, rotation_angle, 'loose', 'bilinear'); % Rotates Image (from calibration script)
  img = imcrop(img, crop_box); % crop according to crop box)
  img = imadjust(img, stretchlim(img),[]); % adjust contrast
  imshow(img); % show image
  
  % Filters  BW threshold

  figure; imshow(img); % show current state of image

  img1 = ~imbinarize(img,black_white_threshold); % make black and white according to threshold
  
  figure; imshow(img1); % show current state of image
  
  % Filters Cleaing the image
  img=img1; clear img1; close all % save image adjustments
  
  se=strel('square',10); % make structural element to filter image pixel size 10
  img1 = imerode(img,se); % removes smaller element around bigger element

  figure; imshow(img1); % show current state of image
  
  close all

  img1 = imfill(img1,'holes'); % fill in holes in large patches 
  se=strel('square',5); % make a finer structural element to filter with pixel size 5
  img1 = imerode(img1,se); % remove smalle elements
  
  % Filters Selecing the sliders
  img1 = bwareafilt(img1, number_of_sliders,'largest'); % Select the two largest objects

  figure; imshow(img1); % show current state of image
  
  % Finding the position of sliders
  
    cen = regionprops(img1, 'centroid'); % Obtains the coordinates to the two patches in the image
  
  % Showing the result
    close all    
    imshow(img1,'InitialMagnification',40); hold on; % show image 
    plot(round(cen(1).Centroid(1)), round(cen(1).Centroid(2)),'rx') %% plot coordinates block 1
    plot(round(cen(2).Centroid(1)), round(cen(2).Centroid(2)),'bx') %% plot coordinates fg
    hold off

%% Section 3 Processing the experiment
close all % close all images
t=[0:2:(length(d)-1)*2]; % make array of time (2 seconds between each image)

n=0; % counter
for time_step=1:number_of_time_steps % for loop over the entire experiment
  n=time_step; % easier to work with n
  
  disp(['Processing: ' d(time_step).name ' ...'])% display image currently worked on in command window
  img = imread([image_folder, filesep, d(time_step).name]); % Load Image
  
  img = img(:,:,2); % Choose one channel of color
  
  img = imrotate(img, rotation_angle, 'loose', 'bilinear'); % Rotates Image (from calibration script)
  img = imcrop(img, crop_box); % crop according to crop box)
  img = imadjust(img, stretchlim(img),[]); % adjust contrast
 % img = imadjust(img, stretchlim(img),[]); % second contrast adjustments

  img = ~imbinarize(img,black_white_threshold); % Convert to black and white 
  
  se=strel('square',10); % create element for filtering
  img = imerode(img,se); % filtering of image
  img = imfill(img,'holes'); % fill holes in aimge
  se=strel('square',5); % create second element for filtering
  img = imerode(img,se); % second filtering of image
  
  img = bwareafilt(img, number_of_sliders,'largest'); % filter image for two largest elements

  cen = regionprops(img, 'centroid'); % get positions of two largest element
  
     for slider = 1:number_of_sliders % save positions of two largest elements in array "p"
        p(n, slider) = fit_x(cen(slider).Centroid(1));           
     end  
  
     % This if statement can be used if one element leaves the image to
     % correctly store coordinates of elements 
%   if  n==1 || (fit_x(cen(number_of_sliders).Centroid(1))-(p(1,number_of_sliders)))>p(n-1,number_of_sliders)
%         disp('save loop 1') 
%       for slider = 1:number_of_sliders
%         p(n, slider) = fit_x(cen(slider).Centroid(1));           
%       end         
%   elseif   (fit_x(cen(number_of_sliders).Centroid(1))-p(1,number_of_sliders))<p(n-1,number_of_sliders)
%        disp('save loop 2'); disp(num2str(n))
%          for slider = 1:number_of_sliders-1
%         p(n,slider+1) = fit_x(cen(slider).Centroid(1));
%         p(n,1)=NaN;
%       end    
%   end
  
if n==1 || mod(n,5)==0 % this if loop shows the image based on a multiple of 5 in order to not sow all images
close all    
imshow(img,'InitialMagnification',40); hold on; 
plot(round(cen(1).Centroid(1)), round(cen(1).Centroid(2)),'rx') % block 1
plot(round(cen(2).Centroid(1)), round(cen(2).Centroid(2)),'bx') % fg / block2 (if three elements)
%plot(round(cen(3).Centroid(1)), round(cen(3).Centroid(2)),'gx')% fg (if three elements)

drawnow % forces the for loop to break and make the image
%pause(1) % can be used if you want the break to be longer to view the
%image
hold off
else  
end

end
disp('Processing Complete')

%% Section 4 Plot
% This sections plots the result for you to check the processing

close all
hold on
h1=plot(t,p(:,1)); % Slider
h2=plot(t,p(:,2)); % Force Gayge
grid on
xlabel('t (s)')
ylabel('x (cm)')
title('Position vs Time')
legend([h1,h2],'Slider','Force gauge','location','SouthEast')
axis tight
box on
hold off


%% Section 5 Preparation and saving
filename=['T_and_P']; % Filename of matlab variable
save(filename,'t','p'); % save output to matlab file


