%% Section 1 Analysis of Stick Slider Data
clear all %remove data
load('T_and_P') % load output from velocimetry

%adjusts position of slider to start from 0
if p(:,1)<0
p(:,1)=p(:,1)+abs(p(1,1)); % block
else
    p(:,1)=p(:,1)-(p(1,1));
end

%adjusts position of force gauge to start from 0
if p(:,2)<0
p(:,2)=p(:,2)+abs(p(1,2));
else
p(:,2)=p(:,2)-(p(1,2)); % fg
end

%% Section 2 Velocity / Loading Rate (from position of force gauge)
close all
t=t'; % transpose t vector

pv=p(:,2); % select position of force gauge
%ind=~isnan(pv); % finds values that arent NaN if needed
v=gradient(pv,t); % calculate gradient i.e. speed

aveV=mean(v); % average speed

%Plot it
plot(v,'k'); hold on; plot([0,length(v)],[aveV,aveV],'-r');
ylim([-1,2])
legend('show','Velocity',[num2str(aveV),' cm/s'],'Location','South')

%% Section 3 Slip v Time Plot

figure % Start figure
set(0,'defaulttextinterpreter','latex') % Set interpreter for plotting
hold on
h1=plot(t,p(:,1),'r-','LineWidth',1); %B1
%h2=plot(t,p(:,2),'b-','LineWidth',1); %fg
%h3=plot(t,p(:,3)','g-','LineWidth',2); % depending
grid on
box on

legend([h1],'b1','fg','Location','SouthEast')
% legend([h1,h2,h3],'fg','b1','b2','Location','SouthEast') % use for more
% blocks

xlim([0 t(end)]) % limit plot from zero to end of time
ylim([0 p(end,2)]) % limit plot from zero to end of positioni

xlabel('t(s)')
ylabel('x(cm)')

title('Slip v Time')
filename=['SlipTime'];
print('-dpdf','-r300',filename) % save to pdf

%% Section 4 Slip length part 1

close all
set(0,'defaulttextinterpreter','latex')
hold on

% take the difference of P array to show sliplengths 
h1=plot(t(1:end-1),diff(p(:,1)),'r-x'); %B1
%h2=plot(t(1:end-1),diff(p(:,2)),'b-','LineWidth',1.5); %fg
%h3=plot(t(1:end-1),diff(p(:,3)),'g-','LineWidth',1.5); %B2 
grid on
box on
legend([h1],'b1','Location','NorthWest')
xlim([0 t(end)]) 
% ylim([0 p(end,1)])

xlabel('t(s)')
ylabel('x(cm)')

title('Slip Lengths v Time')
filename=['SlipLengths'];
print('-dpdf','-r300',filename)


%% Section 4 Slip length part 2 Set threshold for Slip Length
close all
sl=diff(p(:,1));% make difference to find places of slips
sl_fil=sl > 0.25; % set threshold for values which show slip

% Show selection
plot(t(1:end-1),sl,'-x'); grid on; axis tight; hold on
plot(t(sl_fil),sl(sl_fil),'rx')

%% Section 5 Calculating total slip lengths
sl(sl<0.25)=0; % set values below threshold value to be zero

clear sl_a
for i=1:1:length(sl) % Collect slip lengths
   if sl_fil(i)==1      
       sl1=sl(i)-sl(i-1);      
       sl_a(i,1)=sl1;
       sl_a(i,2)=t(i);       
       clear sl1      
   else
   end   
end
clear sl_b
for i=1:1:length(sl) % Sum slip lengths for a given slip event to get total slip length for that event
    % THIS ASSUMUES ONLY FOUR SUBSEQUENT SLIPS
    if i==1
    % do noting
    elseif length(sl)-i<4
    % do nothing    
    elseif sl_a(i-1,1)==0 && sl_a(i,1) ~= 0 && sl_a(i+1,1) == 0  % 1 event    
        sl_b(i,1)=sl_a(i,1);
        sl_b(i,2)=t(i);              
    elseif sl_a(i-1,1)==0 && sl_a(i,1) ~= 0 && sl_a(i+1,1) ~= 0 && sl_a(i+2,1) ==0  % 2 events   
        sl_b(i,1)= sl_a(i,1)+sl_a(i+1,1);
        sl_b(i,2)=t(i);      
   elseif sl_a(i-1,1)==0 && sl_a(i,1) ~= 0 && sl_a(i+1,1) ~= 0 && sl_a(i+2,1) ~=0 % 3 events
       sl_b(i,1)= sl_a(i,1)+sl_a(i+1,1)+sl_a(i+2,1);
       sl_b(i,2)=t(i);
   elseif sl_a(i-1,1)==0 && sl_a(i,1) ~= 0 && sl_a(i+1,1) ~= 0 && sl_a(i+2,1) ~=0 && sl_a(i+3,1) ~=0 % 4 events
       sl_b(i,1)= sl_a(i,1)+sl_a(i+1,1)+sl_a(i+2,1)+sl_a(i+3,1);
       sl_b(i,2)=t(i);           
    else
    end
end
clear sd_a


%% Section 6 Mean, median and standard deviation
close all

sl_f=sl_b(:,1); % Select slip lengths
sl_f(sl_f == 0)=NaN; % Make zeros to NaNs (easier to remove)

ave_sd=nanmean(abs(sl_f)); % mean
med_sd=nanmedian(abs(sl_f)); % median
std_sd=nanstd(abs(sl_f)); % standard deviation
max_sd=max(abs(sl_f)); % max value
min_sd=min(abs(sl_f)); % min value

%histogram plot
hist(abs(sl_f),20); grid on; hold on
plot([ave_sd ave_sd], [0 25],'r');
plot([med_sd med_sd], [0 25],'g');
plot([ave_sd+std_sd ave_sd+std_sd], [0 25],'y');
plot([ave_sd-std_sd ave_sd-std_sd], [0 25],'y');

    
%% Section 7 part 1 Control threshold for what is defined as slip events
s1=diff(p(:,1)); % slip events
L1=s1>0.25; % Slip events thresholded IMPORTANT

%s2=diff(p(:,3)); % For second slider
%L2=s2>0.8;

I3=zeros(size(L1)); 
I3(2:end)=L1(2:end)-L1(1:end-1);
I3(I3<0)=0;
L1=I3; clear I3;

%I3=zeros(size(L2));
%I3(2:end)=L2(2:end)-L2(1:end-1)
%I3(I3<0)=0;
%L2=I3; clear I3;

%logi=L1+L2; % multiple events
%logi=logi>1;

close all
set(0,'defaulttextinterpreter','latex')
hold on

%subplot(3,1,1)
h2=bar(t(1:end-1),cumsum(L1));%B1
xlim([0 t(end)]) 
grid on
box on
ylabel('Slip Event')
legend('b1','Location','NorthWest') 

% subplot(3,1,2)
% h3=bar(t(1:end-1),cumsum(L2),'g'); %B2
%  xlim([200 t(end)]) 
% grid on
% box on
% ylabel('Slip Event') 
% legend('b2','Location','NorthWest') 
% 
% subplot(3,1,3)
% h1=bar(t(1:end-1),cumsum(logi),'r','stacked'); %FG
% xlim([200 t(end)])
% grid on
% box on
% ylabel('Slip Event')
% legend('Common','Location','NorthWest') 
% 
% 
% % legend([h2,h3,h1],'b1','b2','common','Location','NorthWest') 
% 
% % ylim([0 p(end,1)*100])
% 
% xlabel('t(s)')

suptitle('Slip Events v Time')
filename=['SlipEvents'];


%% Section 7 part 1 display slip events

close all
set(0,'defaulttextinterpreter','latex')
hold on

h1=plot(t,p(:,1),'r-','LineWidth',2); %FG
%h2=plot(t,p(:,2)*100,'b-','LineWidth',2); %B
%h3=plot(t,p(:,3)*100,'g-','LineWidth',2); %B2
bar(t(2:end),L1*500,'k')
%bar(t(2:end),L2*1000,'g')
%bar(t(2:end),logi*1000,'r')


grid on
box on
% legend([h1,h2,h3],'fg','b1','b2','Location','SouthEast')
%legend([h2,h3],'b1','b2','Location','SouthEast')

xlim([0 t(end)]) 
ylim([0 p(end,2)])

xlabel('t(s)')
ylabel('x(cm)')


title('Slip v Time & Events')
filename=['SlipTimeEvents'];
print('-dpdf','-r300',filename)
print('-dpdf','-r300',filename)


%% Section 8 Sticking, sliping or creeping? aka velocity of block
close all

pv=p(:,1); % get position and time from p array
ind=~isnan(pv); % find any nans if any

v=gradient(pv,t); % take the gradient i.e. velocity of the block
aveV=mean(v);

plot(t,v,'k'); hold on; axis tight; grid on; 
plot([0,t(end)],[aveV,aveV],'-r');

xlabel('t (s)')
ylabel('x (cm/s)')

%Percentages of behaviour
stick=(sum(v < 0.05)*2)/t(end);
slide=(sum(v > 0.05 & v <0.2)*2)/t(end);
slip=(sum(v>0.2))*2/t(end);

