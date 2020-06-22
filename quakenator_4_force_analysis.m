%% Section 1 Load and Format Data

clear all; clc; close all

% Load data
F=readtable('Force7.txt'); % Read data
F=F.Var1'; % select column 1
t_f=((0:1:numel(F)-1)*0.25); % make array with time

nanf=isnan(F); % Find NaN

plot(t_f,F,'b'); grid on; axis tight; hold on % Plot t and f
plot(t_f(nanf), nanf*nanmean(F),'xr','MarkerSize',10); % plot location of Nan
%% Section 2 interpolation over nan

X = ~isnan(F); % (find nans)
Y = cumsum(X-diff([1,X])/2); % create t to interpolate over
Z = interp1(1:nnz(X),F(X),Y); % inerpolate and nnz is the number of nonzero elements

%% Section 3 Plot to show interpolation results
close all
plot(t_f,Z,'b'); grid on; axis tight; hold on; % show t and f now called Z
plot(t_f(nanf),Z(nanf),'r.') % show where the interpolation is


%% Section 4 force drop analysis starts
close all
sd=diff(Z);% make difference to find places of stress drops
sd_fil=sd < -0.1; % set threshold for values which show stress drops

% Show selection
plot(t_f(1:end-1),sd,'-x'); grid on; axis tight; hold on
plot(t_f(sd_fil),sd(sd_fil),'rx')

sd_fil=logical([0,sd_fil]); % add extra zero at start to make time consistent

close all % Plot to show
plot(t_f,Z,'-x'); grid on; axis tight; hold on
plot(t_f(sd_fil),Z(sd_fil),'x'); grid on; axis tight; 

%% Section 5 calculating total force drops
clear sd_a
for i=1:1:length(Z) % Collect stress drops
   if sd_fil(i)==1      
       sd1=Z(i)-Z(i-1);      
       sd_a(i,1)=sd1;
       sd_a(i,2)=t_f(i);       
       clear sd1      
   else
   end   
end
clear sd_b
for i=1:1:length(sd_a) % summ stress drops for a given slip event to get total stress drop
    % THIS ASSUMUES ONLY THREE SUBSEQUENT DROPS
    if i==1
    % do noting
    elseif length(sd_a)-i<3
    % do nothing    
    elseif sd_a(i-1,1)==0 && sd_a(i,1) ~= 0 && sd_a(i+1,1) == 0  % 1 event    
        sd_b(i,1)=sd_a(i,1);
        sd_b(i,2)=t_f(i);              
    elseif sd_a(i-1,1)==0 && sd_a(i,1) ~= 0 && sd_a(i+1,1) ~= 0 && sd_a(i+2,1) ==0  % 2 events   
        sd_b(i,1)= sd_a(i,1)+sd_a(i+1,1);
        sd_b(i,2)=t_f(i);      
   elseif sd_a(i-1,1)==0 && sd_a(i,1) ~= 0 && sd_a(i+1,1) ~= 0 && sd_a(i+2,1) ~=0 % 3 events
       sd_b(i,1)= sd_a(i,1)+sd_a(i+1,1)+sd_a(i+2,1);
       sd_b(i,2)=t_f(i);     
    else
    end
end
clear sd_a

%% Section 6 Mean, median and standard deviation
close all

sd_f=sd_b(:,1); % select final force drop
sd_f(sd_f == 0)=NaN; % Make all zeros equal NaN

ave_sd=nanmean(abs(sd_f)); % calculate mean
med_sd=nanmedian(abs(sd_f)); % calculate median
std_sd=nanstd(abs(sd_f)); % calculated standard deviation
max_sd=max(abs(sd_f)); % calculate max
min_sd=min(abs(sd_f)); % calculate min

% plot histogram and stats
hist(abs(sd_f),20); grid on; hold on
plot([ave_sd ave_sd], [0 50],'r');
plot([med_sd med_sd], [0 50],'g');
plot([ave_sd+std_sd ave_sd+std_sd], [0 50],'y');
plot([ave_sd-std_sd ave_sd-std_sd], [0 50],'y');

