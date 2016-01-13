%% Combining all data sets
% loading the combined data and making the master matrix
% Currently one has to load them manually from the folder where one has
% saved the values matirx as all of them have the same name one has to do
% it one at a time by deleting the older values matirx and saving these in
% separate variables which can then be concatenated. I haven't yet figured
% out an efficent and less time consuming way than this
load('values.mat');
x17=x(:,:);
y17=y(:,:);
distance17=distance(:,:);
angle17=angle(:,:);
load('1408190953_017_h_3.mat');
Spikes17=V1408190953_017_nw_2.values;
load('1408190953_017_i_3.mat');
Spikes17B=V1408190953_017_nw_2.values;
load('values.mat');
x18=x(:,:);
y18=y(:,:);
distance18=distance(:,:);
angle18=angle(:,:);
load('1408191054_018_h_3.mat');
Spikes18=V1408191054_018_nw_2.values;
load('1408191054_018_i_3.mat');
Spikes18B=V1408191054_018_nw_2.values;
load('values.mat');
x19=x(:,:);
y19=y(:,:);
distance19=distance(:,:);
angle19=angle(:,:);
load('1408191154_019_h_3.mat');
Spikes19=V1408191154_019_nw_2.values;
load('1408191154_019_i_3.mat')
Spikes19B=V1408191154_019_nw_2.values;
load('values.mat');
x20=x(:,:);
y20=y(:,:);
distance20=distance(:,:);
angle20=angle(:,:);
load('1408191255_020_h_3.mat');
Spikes20=V1408191255_020_nw_2.values;
load('1408191255_020_i_3.mat');
Spikes20B=V1408191255_020_nw_2.values;
load('values.mat');
x21=x(:,:);
y21=y(:,:);
distance21=distance(:,:);
angle21=angle(:,:);
load('1408191355_021_h_3.mat');
Spikes21=V1408191355_021_nw_2.values;
load('1408191356_021_i_3.mat');
Spikes21B=V1408191355_021_nw_2.values;
load('values.mat');
x22=x(:,:);
y22=y(:,:);
distance22=distance(:,:);
angle22=angle(:,:);
load('1408191456_022_h_3.mat');
Spikes22=V1408191456_022_Ch3.values;
load('1408191456_022_i_3.mat');
Spikes22B=V1408191456_022_Ch5.values;
load('values.mat');
x23=x(:,:);
y23=y(:,:);
distance23=distance(:,:);
angle23=angle(:,:);
load('1408191557_023_h_3.mat');
Spikes23=V1408191557_023_nw_2.values;
load('1408191557_023_i_3.mat');
Spikes23B=V1408191557_023_nw_2.values;
% Now combine all the data into a bigger matrix.
Spikes=[Spikes17;Spikes18;Spikes19;Spikes20;Spikes21];
Spikes2=[Spikes17B;Spikes18B;Spikes19B;Spikes20B;Spikes21B];
x=[x17(1:36000,:);x18(1:36000,:);x19(1:36000,:);x20(1:36000,:);x21(1:36000,:)];
y=[y17(1:36000,:);y18(1:36000,:);y19(1:36000,:);y20(1:36000,:);y21(1:36000,:)];
distance=[distance17(1:36000,:);distance18(1:36000,:);distance19(1:36000,:);distance20(1:36000,:);distance21(1:36000,:)];
angle=[angle17(1:36000,:);angle18(1:36000,:);angle19(1:36000,:);angle20(1:36000,:);angle21(1:36000,:)];
x=double(x);
y=double(y);
distance=double(distance);
angle=double(angle);
x(100391:100401,1)=NaN;
y(100391:100401,1)=NaN;
%% The first half of the script deals with calculating different parameters of interest the second half of the script is solely for plotting the values calculated here.

% Running rolling average
% For Channel 1
% k=1;
% counter=1;
% spalte=size(Spikes,2)+1;
% while k<size(Spikes,1)
%     
%     if Spikes(k,1)==0
%         
%         k=k+1;
%         if Spikes(k,1)==0
%             counter=counter+1;
%         else
%             Spikes(k-counter:k-1,spalte)=10/counter;
%             counter=1;
%         end
%     else
%         Spikes(k,spalte)=Spikes(k,1)*10;
%         k=k+1;
%     end
% end
% Performs rolling average for Channel1
window=5;
for k=window:size(Spikes,1)-window
    temp=Spikes(k-window+1:k+window,1);
    Spikes(k,2)=10*(sum(temp)./(2*window));
end


% For Channel 2
% k=1;
% counter=1;
% spalte=size(Spikes2,2)+1;
% while k<size(Spikes2,1)
%     
%     if Spikes2(k,1)==0
%         
%         k=k+1;
%         if Spikes2(k,1)==0
%             counter=counter+1;
%         else
%             Spikes2(k-counter:k-1,spalte)=10/counter;
%             counter=1;
%         end
%     else
%         Spikes2(k,spalte)=Spikes2(k,1)*10;
%         k=k+1;
%     end
% end


% Performs rolling average for Channel2
for k=window:size(Spikes2,1)-window
    temp=Spikes2(k-window+1:k+window,1);
    Spikes2(k,2)=10*(sum(temp)./(2*window));
end
% Running the speeds
distance(100391:100401,1)=NaN;
% Speeds=zeros(size(distance,1),12); % Creates an empty matrix to store the variables in the loop
% for i=1:size(distance,2)
%     for j=1:size(distance,1)-1
%         Speeds(j,i)=abs(distance(j+1,i)-distance(j,i))/0.10;
%     end
% end
Speeds=distance;
% For the recorded bee
% Trajectory plot of recorded bee
Recx=x(:,1); % recorded bee x coordinates
Recy=y(:,1);% recorded bee y coordinates
A=Recx;
B=Recy;
A(A==0)=NaN; % Convert Zeros to NaN
B(B==0)=NaN; % Convert Zeros to NaN

% Distribution of orientation angle of recorded bee
Ang=angle(:,1);
Ang(100391:100401)=NaN;
Angr=degtorad(Ang); % Converts the degree to radians.


% Plot of spike frequency per orientation angle
N=Spikes(:,2);
N2=Spikes2(:,2);
Ang_recBee=angle(:,1);% Takes the angle of the recorded bee
binEdge = linspace(min(Ang_recBee),max(Ang_recBee),100); % the bins I want the angles to be divided
[n,bin] = histc(Ang_recBee,binEdge); % n is the sum of all the elements in the bin and bin gives the number where the element belongs to.
% Passing the Spikes over a high and low pass filter respectively.
parameters.RC=300; % RC is the time constant
yy=lowpass(Spikes(:,2),1,parameters); % Passes the Spikes over a low pass filter. For more information check low pass filter function.
yyy=lowpass(Spikes2(:,2),1,parameters);
parameters.RC=30; % RC is the time constant
cc=highpass(Spikes(:,2),1,parameters); % Passes the Spikes over a high pass filter. For more information check high pass filter function.
ccc=highpass(Spikes2(:,2),1,parameters);
delta=cc-ccc; % Difference between the high pass
% Divide the recordings into contacts and then run the entire analysis in blocks
Length_bee= 24.3846; % average lenght of the recorded bee.
SRT=zeros(size(x,1),210); % create the empty matrix to put the different categories based on the distance of the appraching bee to the recorded bee.
RD= zeros(size(x,1),210);
SG= zeros(size(x,1),210);
VVSL= zeros(size(x,1),210);
% The program loops through the script and checks for the 3 conditions of
% the distance apporaching bee with respect to the recorded bee and takes
% into considereation the spike frequencies, coordinate, speed,angle of
% the bee at those time instances where the distance has been divided into
% different classes.
for j=2:12
    for i=1:size(x,1)
        if CalcDistance(x(i,1),y(i,1),x(i,j),y(i,j))<(1.5*Length_bee)
            RD(i,j)=Spikes(i,2);
            RD(i,j+12)=Spikes2(i,2);
            RD(i,j+24)=Speeds(i,j);
            RD(i,j+36)=angle(i,j);
            RD(i,j+48)=x(i,j);
            RD(i,j+60)=y(i,j);
            RD(i,j+72)=angle(i,1)-angle_coordinate(x(i,1),y(i,1),x(i,j),y(i,j));
            RD(i,j+84)=times(i,1);
            RD(i,j+96)=yy(i,1);
            RD(i,j+108)=yyy(i,1);
            RD(i,j+120)=cc(i,1);
            RD(i,j+132)=ccc(i,1);
            RD(i,j+144)=delta(i,1);
        elseif CalcDistance(x(i,1),y(i,1),x(i,j),y(i,j))>(1.5*Length_bee)&& CalcDistance(x(i,1),y(i,1),x(i,j),y(i,j))<(5*Length_bee)
            VVSL(i,j)=Spikes(i,2);
            VVSL(i,j+12)=Spikes2(i,2);
            VVSL(i,j+24)=Speeds(i,j);
            VVSL(i,j+36)=angle(i,j);
            VVSL(i,j+48)=x(i,j);
            VVSL(i,j+60)=y(i,j);
            VVSL(i,j+72)=angle(i,1)-angle_coordinate(x(i,1),y(i,1),x(i,j),y(i,j));
            VVSL(i,j+96)=yy(i,1);
            VVSL(i,j+108)=yyy(i,1);
            VVSL(i,j+120)=cc(i,1);
            VVSL(i,j+132)=ccc(i,1);
            VVSL(i,j+144)=delta(i,1);
            VVSL(i,j+84)=times(i,1);
        elseif CalcDistance(x(i,1),y(i,1),x(i,j),y(i,j))>(5*Length_bee) && CalcDistance(x(i,1),y(i,1),x(i,j),y(i,j))<(10*Length_bee)
            SG(i,j)=Spikes(i,2);
            SG(i,j+12)=Spikes2(i,2);
            SG(i,j+24)=Speeds(i,j);
            SG(i,j+36)=angle(i,j);
            SG(i,j+48)=x(i,j);
            SG(i,j+60)=y(i,j);
            SG(i,j+72)=angle(i,1)-angle_coordinate(x(i,1),y(i,1),x(i,j),y(i,j));
            SG(i,j+84)=times(i,1);
            SG(i,j+96)=yy(i,1);
            SG(i,j+108)=yyy(i,1);
            SG(i,j+120)=cc(i,1);
            SG(i,j+132)=ccc(i,1);
            SG(i,j+144)=delta(i,1);
        else
            SRT(i,j)=Spikes(i,2);
            SRT(i,j+12)=Spikes2(i,2);
            SRT(i,j+24)=Speeds(i,j);
            SRT(i,j+36)=angle(i,j);
            SRT(i,j+48)=x(i,j);
            SRT(i,j+60)=y(i,j);
            SRT(i,j+72)=angle(i,1)-angle_coordinate(x(i,1),y(i,1),x(i,j),y(i,j));
            SRT(i,j+84)=times(i,1);
            SRT(i,j+96)=yy(i,1);
            SRT(i,j+108)=yyy(i,1);
            SRT(i,j+120)=cc(i,1);
            SRT(i,j+132)=ccc(i,1);
            SRT(i,j+144)=delta(i,1);
        end
        
    end
end
% length of interactions
% for this create a counter and count the maximum frames before the
% counter resets to 1, this would give the user, the time length in number of frames for each
% interaction.
AK=zeros(size(x,1),500);
BK=zeros(size(x,1),500);
CK=zeros(size(x,1),500);
DK=zeros(size(x,1),500);
for j=2:12
    for i=2:size(x,1)
        if abs(RD(i-1,j+84)-RD(i,j+84))>1
            AK(i,j)=1;
        else
            AK(i,j)=AK(i-1,j)+1;
        end
    end
end
for j=2:12
    for i=2:size(x,1)
        if abs(SRT(i-1,j+84)-SRT(i,j+84))>1
            DK(i,j)=1;
        else
            DK(i,j)=DK(i-1,j)+1;
        end
    end
end
for j=2:12
    for i=2:size(x,1)
        if abs(SG(i-1,j+84)-SG(i,j+84))>1
            CK(i,j)=1;
        else
            CK(i,j)=CK(i-1,j)+1;
        end
    end
end
for j=2:12
    for i=2:size(x,1)
        if abs(VVSL(i-1,j+84)-VVSL(i,j+84))>1
            BK(i,j)=1;
        else
            BK(i,j)=BK(i-1,j)+1;
        end
    end
end
% Here,the length of each interaction is found along with the mean speed
% and mean angle of orientation of the apporaching bee.
for j=2:12
    Spikes1=RD(:,j);
    Spikes11=RD(:,j+12);
    speed=RD(:,j+24);
    Angles=RD(:,j+72);
    HHH=find(AK(:,j)==1);
    
    for i=1:length(HHH)
        if i==1
            V=AK(:,j);
            AK(i,j+12)=length(V(1:HHH(i)));
            AK(i,j+24)=mean(speed(1:HHH(i)));
            AK(i,j+36)=mean(Angles(1:HHH(i)));
            AK(i,j+48)=mean(Spikes1(1:HHH(i)));
            AK(i,j+60)=mean(Spikes11(1:HHH(i)));
        else
            V=AK(:,j);
            AK(i,j+12)=length(V(HHH(i-1):HHH(i)));
            AK(i,j+24)=mean(speed(HHH(i-1):HHH(i)));
            AK(i,j+36)=mean(Angles(HHH(i-1):HHH(i)));
             AK(i,j+48)=mean(Spikes1(HHH(i-1):HHH(i)));
            AK(i,j+60)=mean(Spikes11(HHH(i-1):HHH(i)));
        end
        
    end
end
for j=2:12
    
    HHHH=find(BK(:,j)==1);
    for i=1:length(HHHH)
        if i==1
            V=BK(:,j);
            BK(i,j+12)=length(V(1:HHHH(i)));
        else
            V=BK(:,j);
            BK(i,j+12)=length(V(HHHH(i-1):HHHH(i)));
        end
        
    end
end
for j=2:12
    HHHHH=find(CK(:,j)==1);
    for i=1:length(HHHHH)
        if i==1
            V=CK(:,j);
            CK(i,j+12)=length(V(1:HHHHH(i)));
        else
            V=CK(:,j);
            CK(i,j+12)=length(V(HHHHH(i-1):HHHHH(i)));
        end
        
    end
end

for j=2:12
    HHHHHH=find(DK(:,j)==1);
    for i=1:length(HHHHHH)
        if i==1
            V=DK(:,j);
            DK(i,j+12)=length(V(1:HHHHHH(i)));
        else
            V=DK(:,j);
            DK(i,j+12)=length(V(HHHHHH(i-1):HHHHHH(i)));
        end
        
    end
end
%% Plotting this part of the script deals with plotting the values calculated in the above scripts.
% Trajectory of recorded bee
% A(100391:100401)=NaN;
% B(100391:100401)=NaN;
figure(1)
plot(A,B)% Plots the X-Y coordinates of the recorded bee
xlim([0 1600])
ylim([0 1200])
xlabel('Xcoordinate of Recorded Bee');
ylabel('Ycoordinate of Recorded Bee');
title('Trajectory of the recorded Bee');
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Trajectory.jpg')
        close all
    case 'NO'
        close all
        
end
% Distribution of orientation angle of recorded bee
figure(2)
h = rose(Angr,200); % rose actually plots a polar plot, and the 20 is the bins of the angles
% x = get(h,'Xdata');
% y = get(h,'Ydata');
% g = patch(x,y,'r');
title('Distribution of orientation angle');

choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('orientation angle recorded bee distribution.jpg')
        close all
    case 'NO'
        
        close all
        
end

% Plotting spikes as a false color on the trajectory.
figure(3)
color_line3(A, B, 1:size(x,1),Spikes(:,2))
% what this function does is it uses suface plots and plots the x and y coordinate with respect
% to the third coordinate which is time and then uses the values of spike
% frequency and color codes the graph.
xlabel('X coordinate of recorded Bee');
ylabel('Y coordinate of recorded Bee');
c = colorbar;
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('fequency_color_trajectory channel 1.jpg')
        close all
    case 'NO'
        close all
        
end

figure(4)
color_line3(A, B, 1:size(x,1),Spikes2(:,2))
% what this function does is it uses suface plots and plots the x and y coordinate with respect
% to the third coordinate which is time and then uses the values of spike
% frequency and color codes the graph.
xlabel('X coordinate of recorded Bee');
ylabel('Y coordinate of recorded Bee');
c = colorbar;
caxis([0 110])

%set(get(c,'title'),'string','Spike Frequency of Channel 1');
%this creates the lable and sets it horizontally on the colour bar
ylabel(c, 'Spike Frequency of Channel 1') % This places it neatly at the side of the colour bar.
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('fequency_color_trajectory channel2.jpg')
        close all
    case 'NO'
        close all
        
end
% Distribution of spike frequency per angle bin
figure(5)
subplot(2,1,1)
boxplot(N2,bin);
h=findobj(gca,'tag','Outliers');
delete(h);
xlabel('groups or bins');
ylabel('spike frequency in HZ');
title('Distribution of spike frequency per angle bin')
subplot(2,1,2)
boxplot(N,bin)
h=findobj(gca,'tag','Outliers');
delete(h);
xlabel('groups or bins');
ylabel('spike frequency in HZ');
title('Distribution of spike frequency per angle bin')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Distibution_Spikes_per_angle_bin.jpg')
        close all
    case 'NO'
        close all
end
% plot of channel 1 vs channel 2
figure(999)
scatter(1:size(x,1),Spikes(:,2))
hold on
scatter(1:size(x,1),Spikes2(:,2),'r')
figure(6)
scatter(Spikes(:,2),Spikes2(:,2),'.')
xlabel('Channel1')
ylabel('Channel2')
title('Plot of Channel 1 vs Channel 2')
axis equal
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Channel1vsChannel2_roll_average.jpg')
        close all
    case 'NO'
        close all
end
% Distribution of Spike frequency of Channel 1 and Channel 2
figure(7)
subplot(2,1,1)
hist(Spikes(:,2),1000)
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of spike frequency of Channel 1')
subplot(2,1,2)
hist(Spikes2(:,2),1000)
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of spike frequency of Channel 2')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Distibution_Spikes_Both_Channels.jpg')
        close all
    case 'NO'
        close all
end
% Plotting the high passed and low passed spikes over the original
figure(8)
plot(1:size(x,1),Spikes(:,2))
hold on
plot(1:size(x,1),yy,'r')
legend('original data','lowpassed data unit 1')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('lowpassed_Spike1.jpg')
        close all
    case 'NO'
        close all
end
figure(9)
plot(1:size(x,1),Spikes2(:,2))
hold on
plot(1:size(x,1),yyy,'r')
legend('original data','lowpassed data unit 2')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('lowpassed_Spike2')
        close all
    case 'NO'
        close all
end
figure(10)
plot(1:size(x,1),Spikes(:,2))
hold on
plot(1:size(x,1),cc,'r')
legend('original data','high passed data of unit 1')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Highpassed_Spike1.jpg')
        close all
    case 'NO'
        close all
end

figure(11)
plot(1:size(x,1),Spikes2(:,2))
hold on
plot(1:size(x,1),ccc,'r')
legend('original data','high passed data of unit 2')
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Highpasses_Spike2.jpg')
        close all
    case 'NO'
        close all
end
% Plotting the both the spikes and speed as the false color for different
% contacts
figure(12)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(RD(:,j), RD(:,j+12), 10, RD(:,j+24),'.');
    colorbar
    title('touch')
end
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Spikes_Spikes2_speed_touch.jpg')
        close all
    case 'NO'
        close all
end
figure(13)
for j=1:12
    %cmp = jet(T); % create the color maps changed as in jet color map
    scatter(SRT(:,j), SRT(:,j+12), 10, SRT(:,j+24), '.');
    colorbar
    title('far')
end
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Spike1_Spike2_Speed_far.jpg')
        close all
    case 'NO'
        close all
end
figure(14)
for j=1:12
    %cmp = jet(L); % create the color maps changed as in jet color map
    scatter(VVSL(:,j), VVSL(:,j+12), 10, VVSL(:,j+24),'filled');
    colorbar
    title('very near')
end
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Spike1_Spike2_Very_near.jpg')
        close all
    case 'NO'
        close all
end
figure(15)
for j=1:12
    %cmp = jet(G); % create the color maps changed as in jet color map
    scatter(SG(:,j), SG(:,j+12), 10, SG(:,j+24), 'filled');
    colorbar
    title('near')
end
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Spike1_Spike2_Speed_near.jpg')
        close all
    case 'NO'
        close all
end
figure(333)
for j=1:12
    hold all
    scatter(SRT(:,j), SRT(:,j+12), '.','r');
    scatter(RD(:,j), RD(:,j+12),'.','g');
    scatter(VVSL(:,j), VVSL(:,j+12) ,'.','b');
    scatter(SG(:,j), SG(:,j+12), '.','y');
    xlabel('spikes of channel 1')
    ylabel('spikes of channel 2')
    legend('far','touch','verynear','near')
end
% Plotting the Distribution of angle of approach of near bees to recorded bee.
figure(16)
hold all
hist(RD(:,73),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,74),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,75),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,76),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,77),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');


hist(RD(:,78),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');



hist(RD(:,79),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,80),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,81),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,82),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');

hist(RD(:,83),20)
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
choice = questdlg('Would you like to save?',...
    'Save Menu',...
    'Yes','NO','NO');
% Handle response
switch choice
    case 'Yes'
        savefig('Distribution of angle of approach of near bees.jpg')
        close all
    case 'NO'
        close all
end


% Distribution of time
V1=AK(:,14:24);
V2=BK(:,14:24);
V3=CK(:,14:24);
V4=DK(:,14:24);
V1(V1==0)=NaN;
V2(V2==0)=NaN;
V3(V3==0)=NaN;
V4(V4==0)=NaN;
figure(17)
hist(V1,50000)
xlim([0 35])
figure(18)
hist(V2)
figure(19)
hist(V3)
figure(20)
hist(V4)
% Spike for closest bee
hist(Y);
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
hold on
hist(PB)
set(get(gca,'child'),'FaceColor','none','EdgeColor','b');
hist(PB,Y)
subplot(2,1,1)
scatter(Y,PBf,'.');
subplot(2,1,2)
scatter(Y,PB,'.')


scatter(Y,PBf,'r','.');
hold on
sg=scatter(Y,PB,'.','g');
sg=patch(Y,PBf,'r','FaceColor','none');
alpha(sg,.1);
% sMarkers=sg.MarkerHandle; %hidden marker handle
% sMarkers.FaceColorData = uint8(255*[1;0;0;0.1]); %fourth element allows setting alpha
% sMarkers.EdgeColorData = uint8(255*[1;0;0;0]); %set edge color in a similar way

% plotting with delta(The difference between the high pass
% filters),highpass,lowpass filters
% false color on trajectory
figure(22)
color_line3(A, B, 1:size(x,1),delta);
c=colorbar;
caxis([-10 20])
xlabel('x coordinate of recorded bee');
ylabel('y coordinate of recorded bee');
title('x-y coordinate with the difference of high pass as a false color')
figure(23)
color_line3(A, B, 1:size(x,1),cc);
c=colorbar;
caxis([-10 20])
xlabel('x coordinate of recorded bee');
ylabel('y coordinate of recorded bee');
title('x-y coordinate with the highpass of unit 1 as a false color')
figure(24)
color_line3(A, B, 1:size(x,1),ccc);
c=colorbar;
caxis([-10 20])
xlabel('x coordinate of recorded bee');
ylabel('y coordinate of recorded bee');
title('x-y coordinate with high pass of unit 2 as a false color')
figure(25)
color_line3(A, B, 1:size(x,1),yy);
c=colorbar;
caxis([-10 20])
xlabel('x coordinate of recorded bee');
ylabel('y coordinate of recorded bee');
title('x-y coordinate with the lowpass of unit 1 as a false color')
figure(26)
color_line3(A, B, 1:size(x,1),yyy);
c=colorbar;
caxis([-10 20])
xlabel('x coordinate of recorded bee');
ylabel('y coordinate of recorded bee');
title('x-y coordinate with the lowpass of unit 2 as a false color')
% Distribution of spike frequency of high pass per bin
figure(27)
subplot(2,1,1)
boxplot(cc,bin);
h=findobj(gca,'tag','Outliers');
delete(h);
xlabel('groups or bins');
ylabel('spike frequency in HZ');
title('Distribution of spike frequency of high pass unit  per angle bin')
subplot(2,1,2)
boxplot(ccc,bin)
h=findobj(gca,'tag','Outliers');
delete(h);
figure(28)
subplot(2,1,1)
boxplot(yy,bin);
h=findobj(gca,'tag','Outliers');
delete(h);
xlabel('groups or bins');
ylabel('spike frequency in HZ');
title('Distribution of spike frequency of low pass per angle bin')
subplot(2,1,2)
boxplot(yyy,bin)
h=findobj(gca,'tag','Outliers');
delete(h);
figure(29)
boxplot(delta,bin);
h=findobj(gca,'tag','Outliers');
delete(h);
xlabel('groups or bins');
ylabel('spike frequency in HZ');
title('Distribution of difference of spike frequency of high pass per angle bin')
%Distribution of Spike frequency of Channel 1 and Channel 2
figure(30)
subplot(2,1,1)
hist(yy,500)
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of low pass spike frequency of unit1')
subplot(2,1,2)
hist(yyy,500)
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of low pass spike frequency of unit2')
figure(31)
subplot(2,1,1)
hist(cc,500)
xlim([-20 20])
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of high pass spike frequency of unit1')
subplot(2,1,2)
hist(ccc,500)
xlim([-20 20])
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of high pass spike frequency of unit2')
figure(32)
hist(delta,500)
set(get(gca,'child'),'FaceColor','none','EdgeColor','b');
xlim([-20 20])
xlabel('frequency in HZ')
ylabel('number of times')
title('distribution of the difference of high pass spike frequency')
% ND=[-20:1/900:20];
% hold on
% norm=normpdf(ND,0,10);
% plot(ND,norm,'r')
% Plotting the both the spikes and speed as the false color for different
% contacts
figure(33)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(RD(:,j+96), RD(:,j+108), 10, RD(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('touch low pass')
end
figure(34)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SRT(:,j+96), SRT(:,j+108), 10, 1:size(SRT,1),'.');
    colorbar
    caxis([1 100])
    title('far low pass')
end
figure(35)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(VVSL(:,j+96), VVSL(:,j+108), 10, VVSL(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('very near low pass')
end
figure(36)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SG(:,j+96), SG(:,j+108), 10, SG(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('near low pass')
end
figure(37)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(RD(:,j+120), RD(:,j+132), 10, RD(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('touch high pass')
end
figure(38)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SRT(:,j+120), SRT(:,j+132), 10, SRT(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('far high pass')
end
figure(39)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(VVSL(:,j+120), VVSL(:,j+132), 10, VVSL(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('very near high pass')
end
figure(40)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SG(:,j+120), SG(:,j+132), 10, SG(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('near high pass')
end
figure(41)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(RD(:,j+84), RD(:,j+144), 10, RD(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('touch delta')
end
figure(42)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SRT(:,j+84), SRT(:,j+144), 10, SRT(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('far delta')
end
figure(43)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(VVSL(:,j+84), VVSL(:,j+144), 10, VVSL(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('very near delta')
end
figure(44)
for j=1:12
    %     cmp = jet(D); % create the color maps changed as in jet color map
    scatter(SG(:,j+84), SG(:,j+144), 10, SG(:,j+84),'.');
    colorbar
    caxis([1 100])
    title('near delta')
end
figure(99)
color_line3(A,B,1:180000,Speeds(:,1))
caxis([0 50])
%% PCA 
Duration=zeros(size(x,1),1); % Create an empty matrix to store all the values of time, speed and angle of al the 12 bees as single column matrix
SP=zeros(size(x,1),1);
AI= zeros(size(x,1),1);
S1=zeros(size(x,1),1);
S2=zeros(size(x,1),1);
for j=2:12
    for i=1:length(RD)
        if j==2
            Duration(i,1)=AK(i,j+12);
            SP(i,1)=AK(i,j+24);
            AI(i,1)=AK(i,j+36);
            S1(i,1)=AK(i,j+48);
            S2(i,1)=AK(i,j+60);
        else
            Duration(i+(size(AK,1)*(j-2)),1)=AK(i,j+12);
            SP(i+(size(AK,1)*(j-2)),1)=AK(i,j+24);
            AI(i+(size(AK,1)*(j-2)),1)=AK(i,j+36);
            S1(i+(size(AK,1)*(j-2)),1)=AK(i,j+48);
            S2(i+(size(AK,1)*(j-2)),1)=AK(i,j+60);
        end
    end
end
MASTERMATRIX=[Duration,SP,AI]; % creates a bigger matix with all the categories
% as different columns
MASTERMATRIX(MASTERMATRIX==0)=NaN;
scatter3(MASTERMATRIX(:,1),MASTERMATRIX(:,2),MASTERMATRIX(:,3))
figure(999)
S1(S1==0)=NaN;
 scatter3(MASTERMATRIX(:,1),MASTERMATRIX(:,2),MASTERMATRIX(:,3),10,S1(:,1));
 xlim([0 1000])
 ylim([0 40])
 %zlim([0 10])
