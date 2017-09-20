% This script is written by Chao-Yuan Yeh. All copyrights are rerserved.
%% Step 1: Choose center reference point, choroid boundary and auxiliary line
clear
close all
% Set scaling factor for OCT image. The scaling factor is specific to every
% OCT machine and should be acquired from the manufacturer.
x_scale = 8/1019*1000;
y_scale = 2.375/768*1000;

[fileName, pathName] = uigetfile('*.tif','Choose OCT image to process');
cd(pathName);
mkdir('MATLAB');
image = imread(fullfile(pathName,fileName));
% create extra margin for the image 
template = ones(868,1219) * double(min(image(:))+120);
template(51:818,101:1119) = image;
image = template; 
clear template
hfig = figure; 
imshow(image,[]);
hfig.Position = [10 10 flip(size(image),2)];

% Draw reference curve for retina
dia = dialog('Position',[100 400 250 100],'Name','Instruction');
uicontrol('Parent',dia,'Style','text',...
    'Position',[10 40 100 40],'String','Draw left reference curve');
uicontrol('Parent',dia,'Position',[85 20 70 25],...
               'String','OK','Callback','delete(gcf)');
clear dia
h = imfreehand('closed',false);wait(h);
linepos = getPosition(h);
delete(h)
clear h;

% Generate a smooth curve from the reference curve drawn by user
xmin = min(linepos(:,1));
xmax = max(linepos(:,1));
x = xmin:0.1:xmax;
p = polyfit(linepos(:,1),linepos(:,2),4);
y = polyval(p,x);
hold on
plot(x,y,'y','linewidth',3)
clear linepos
curveDist = sqrt(((x(2:end)-x(1:end-1))*x_scale).^2 + ((y(2:end)-y(1:end-1))*y_scale).^2);
% It is important to bear in mind that imshow method displays the image
% with y direction pointing downward ( axis ij ). 
slope = diff(y)/0.1;
islope = -ones(size(slope))./slope;
theta = atan(islope);

% Draw choroid contour
dia = dialog('Position',[100 400 250 100],'Name','Instruction');
uicontrol('Parent',dia,'Style','text',...
    'Position',[10 40 100 40],'String','Draw left choroid contour');
uicontrol('Parent',dia,'Position',[85 20 70 25],...
               'String','OK','Callback','delete(gcf)');
clear dlg

h = imfreehand;
wait(h);
Mask = single(createMask(h));
choroidpos = bwboundaries(Mask); % This method will get equally spaced points "choroidpos = getPopsition(h);" won't.
delete(h);
clear Mask h;
choroidpos = choroidpos{1};
choroidpos = flip(choroidpos,2);
base = length(choroidpos);
bq = 1:.25:base;
ROIX = interp1((1:base),choroidpos(:,1),bq,'spline'); 
ROIY = interp1((1:base),choroidpos(:,2),bq,'spline');
choroidpos = cat(2,squeeze(ROIX)',squeeze(ROIY)');

clear base bq ROIX ROIY
%% Step2 Rotation method
% In this method, when determining the choroid thickness sectioned by each
% reference line, the choroid ROI positions were rotated so that the
% reference line is at ([0 0], [0 200]). This will avoid distance
% calculation between every choroid ROI position and every point on the
% reference line. Instead, distance to reference line is simply the
% absolute value of the x coordinate after rotation. Compared with the
% method that calculates distance point-by-point, this one has 300 fold
% increase in speed. 

extension = 200;
sequence = ( ( 1:length(1:10:length(x)-1) ) -1)*10+1;
x_subset = x(sequence);
y_subset = y(sequence);
islope_subset = islope(sequence);
theta_subset = theta(sequence);
clear data

for ii = 1:length(sequence) 
    
    if islope_subset(ii) < 0
        angleR = .5*pi+theta_subset(ii);
        rotationMatrix = [cos(-angleR) -sin(-angleR) ; sin(-angleR) cos(-angleR)];
        rotationMatrix_R = [cos(angleR) -sin(angleR) ; sin(angleR) cos(angleR)];
    elseif islope_subset(ii) > 0 
        angleR = .5*pi-theta_subset(ii);
        rotationMatrix = [cos(angleR) -sin(angleR) ; sin(angleR) cos(angleR)];
        rotationMatrix_R = [cos(-angleR) -sin(-angleR) ; sin(-angleR) cos(-angleR)];
    else
        rotationMatrix = [1 0 ; 0 1];
    end
    
    % Rotate around the origin of each refence line
    choroidpos_r = mtimes(rotationMatrix, (choroidpos +...
        cat(2, ones(length(choroidpos),1)*(-x_subset(ii)), ones(length(choroidpos),1)*(-y_subset(ii))))')';
    comp = choroidpos_r(:,1) > -1 & choroidpos_r(:,1) <1 &  choroidpos_r(:,2) >0 &  choroidpos_r(:,2) < extension;
    if sum(comp(:)) > 2
        sortedSubset = choroidpos_r(comp,:);
        posind1 = sortedSubset(sortedSubset(:,2) > mean(sortedSubset(:,2)) & ...
            abs(sortedSubset(:,1)) == min(abs(sortedSubset(sortedSubset(:,2) > mean(sortedSubset(:,2)),1))),:)';
        posind1 = posind1(:,1);
        posind2 = sortedSubset(sortedSubset(:,2) < mean(sortedSubset(:,2)) & ...
            abs(sortedSubset(:,1)) == min(abs(sortedSubset(sortedSubset(:,2) < mean(sortedSubset(:,2)),1))),:)';
        posind2 = posind2(:,1);
        posind1 = mtimes(rotationMatrix_R,posind1) + [x_subset(ii);y_subset(ii)];
        posind2 = mtimes(rotationMatrix_R,posind2) + [x_subset(ii);y_subset(ii)];
        data(ii).x1 = posind1(1);
        data(ii).y1 = posind1(2);
        data(ii).x2 = posind2(1);
        data(ii).y2 = posind2(2);
        index = sequence(ii);
        data(ii).dist = -sum(curveDist(index:end));
        data(ii).thickness = sqrt(((posind1(1)-posind2(1))*x_scale)^2 + ((posind1(2)-posind2(2))*y_scale)^2);
    end
end

for ii = 1:length(sequence)
    if islope_subset(ii) < 0
        line([x_subset(ii), x_subset(ii) - extension*cos(theta_subset(ii))],...
            [y_subset(ii), y_subset(ii) - extension*sin(theta_subset(ii))], 'color','r');
    elseif islope_subset(ii) > 0 
        line([x_subset(ii), x_subset(ii) + extension*cos(theta_subset(ii))],...
            [y_subset(ii), y_subset(ii) + extension*sin(theta_subset(ii))], 'color','r');
    else
        line([x_subset(ii), x_subset(ii)],...
            [y_subset(ii), y_subset(ii) + extension], 'color','r');
    end
end

data_L = struct2array(data); % This will remove empty entry. 
data_L = reshape(data_L,6,numel(data_L)/6)';
clear data
for ii = 1:length(data_L) 
    line([data_L(ii,1),data_L(ii,3)],[data_L(ii,2),data_L(ii,4)],'color','g')
end

%% Repeat the process for right side
% Draw reference curve and generate fitted points
dia = dialog('Position',[100 400 250 100],'Name','Instruction');
uicontrol('Parent',dia,'Style','text',...
    'Position',[10 40 100 40],'String','Draw right reference curve');
uicontrol('Parent',dia,'Position',[85 20 70 25],...
               'String','OK','Callback','delete(gcf)');
clear dia

h = imfreehand('closed',false);wait(h);
linepos = getPosition(h);
delete(h)
clear h;
%linepos = sortrows(linepos,1);

xmin = min(linepos(:,1));
xmax = max(linepos(:,1));
x = xmin:0.1:xmax;
p = polyfit(linepos(:,1),linepos(:,2),4);
y = polyval(p,x);
plot(x,y,'y','linewidth',3)
clear linepos

curveDist = sqrt(((x(2:end)-x(1:end-1))*x_scale).^2 + ((y(2:end)-y(1:end-1))*y_scale).^2);
slope = diff(y)/0.1;
islope = -ones(size(slope))./slope;
theta = atan(islope);

% Draw choroid contour
dia = dialog('Position',[100 400 250 100],'Name','Instruction');
uicontrol('Parent',dia,'Style','text',...
    'Position',[10 40 100 40],'String','Draw right choroid contour');
uicontrol('Parent',dia,'Position',[85 20 70 25],...
               'String','OK','Callback','delete(gcf)');
clear dlg

h = imfreehand;wait(h);
Mask = single(createMask(h));
choroidpos = bwboundaries(Mask); % This method will get equally spaced points "choroidpos = getPopsition(h);" won't.
delete(h);
clear Mask h;
choroidpos = choroidpos{1};
%choroidpos = flip(choroidpos,2);
choroidpos = flipdim(choroidpos,2);

base = length(choroidpos);
bq = 1:.25:base;
ROIX = interp1((1:base),choroidpos(:,1),bq,'spline'); 
ROIY = interp1((1:base),choroidpos(:,2),bq,'spline');
choroidpos = cat(2,squeeze(ROIX)',squeeze(ROIY)');
clear base bq ROIX ROIY

extension = 200;
sequence = ( ( 1:length(1:10:length(x)-1) ) -1)*10+1;
x_subset = x(sequence);
y_subset = y(sequence);
islope_subset = islope(sequence);
theta_subset = theta(sequence);

tic
for ii = 1:length(sequence) 
    
    if islope_subset(ii) < 0
        angleR = .5*pi+theta_subset(ii);
        rotationMatrix = [cos(-angleR) -sin(-angleR) ; sin(-angleR) cos(-angleR)];
        rotationMatrix_R = [cos(angleR) -sin(angleR) ; sin(angleR) cos(angleR)];
    elseif islope_subset(ii) > 0 
        angleR = .5*pi-theta_subset(ii);
        rotationMatrix = [cos(angleR) -sin(angleR) ; sin(angleR) cos(angleR)];
        rotationMatrix_R = [cos(-angleR) -sin(-angleR) ; sin(-angleR) cos(-angleR)];
    else
        rotationMatrix = [1 0 ; 0 1];
    end
    
    % Rotate around the origin of each refence line
    choroidpos_r = mtimes(rotationMatrix, (choroidpos +...
        cat(2, ones(length(choroidpos),1)*(-x_subset(ii)), ones(length(choroidpos),1)*(-y_subset(ii))))')';
    comp = choroidpos_r(:,1) > -1 & choroidpos_r(:,1) <1 &  choroidpos_r(:,2) >0 &  choroidpos_r(:,2) < extension;
    if sum(comp(:)) > 2
        sortedSubset = choroidpos_r(comp,:);
        posind1 = sortedSubset(sortedSubset(:,2) > mean(sortedSubset(:,2)) & ...
            abs(sortedSubset(:,1)) == min(abs(sortedSubset(sortedSubset(:,2) > mean(sortedSubset(:,2)),1))),:)';
        posind1 = posind1(:,1);
        posind2 = sortedSubset(sortedSubset(:,2) < mean(sortedSubset(:,2)) & ...
            abs(sortedSubset(:,1)) == min(abs(sortedSubset(sortedSubset(:,2) < mean(sortedSubset(:,2)),1))),:)';
        posind2 = posind2(:,1);
        posind1 = mtimes(rotationMatrix_R,posind1) + [x_subset(ii);y_subset(ii)];
        posind2 = mtimes(rotationMatrix_R,posind2) + [x_subset(ii);y_subset(ii)];
        data(ii).x1 = posind1(1);
        data(ii).y1 = posind1(2);
        data(ii).x2 = posind2(1);
        data(ii).y2 = posind2(2);
        index = sequence(ii);
        data(ii).dist = sum(curveDist(1:index));
        data(ii).thickness = sqrt(((posind1(1)-posind2(1))*x_scale)^2 + ((posind1(2)-posind2(2))*y_scale)^2);
    end
end
toc
for ii = 1:length(sequence)
    if islope_subset(ii) < 0
        line([x_subset(ii), x_subset(ii) - extension*cos(theta_subset(ii))],...
            [y_subset(ii), y_subset(ii) - extension*sin(theta_subset(ii))], 'color','r');
    elseif islope_subset(ii) > 0 
        line([x_subset(ii), x_subset(ii) + extension*cos(theta_subset(ii))],...
            [y_subset(ii), y_subset(ii) + extension*sin(theta_subset(ii))], 'color','r');
    else
        line([x_subset(ii), x_subset(ii)],...
            [y_subset(ii), y_subset(ii) + extension], 'color','r');
    end
end
data_R = struct2array(data); % This will remove empty entry. 
data_R = reshape(data_R,6,numel(data_R)/6)';
clear data

for ii = 1:length(data_R) 
    line([data_R(ii,1),data_R(ii,3)],[data_R(ii,2),data_R(ii,4)],'color','g')
end

%% Post-processing data
data = cat(1, data_L, data_R);
data(:,5) = data(:,5)/8; % scale distance so interval is close to 1.
temp = data(:,5:6);
temp(:,1) = floor(temp(:,1));
temp2 = diff(temp(:,1));
temp(abs(temp2)==0,:) = [];
temp2 = diff(temp(:,1));
index = find(abs(temp2)==2);
index2 = index + (1:length(index))';
temp3 = zeros(size(temp,1)+length(index),size(temp,2));
temp3(index2,:) = nan;
temp3(~isnan(temp3))=temp;
temp3(index2,:) = (temp3(index2-1,:) + temp3(index2+1,:))/2;
disData = temp3;
figure; scatter(disData(:,1), disData(:,2))

%% Save processed data
disData(:,3) = disData(:,2);
disData(:,2) = 0;

cd('MATLAB')
if strfind(fileName,'OD_DC_0_Degree_1')
    data_OD_DC_0_Degree_1 = disData;
    save('OD_DC_0_Degree_1.mat','data_OD_DC_0_Degree_1');
elseif strfind(fileName,'OD_DC_0_Degree_2')
    rotationMatrix = [cos(pi/2) -sin(pi/2) ; sin(pi/2) cos(pi/2)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    data_OD_DC_0_Degree_2 = disData;
    save('OD_DC_0_Degree_2.mat','data_OD_DC_0_Degree_2');
elseif strfind(fileName,'OD_DC_45_Degree_1')
    rotationMatrix = [cos(-pi/4) -sin(-pi/4) ; sin(-pi/4) cos(-pi/4)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    data_OD_DC_45_Degree_1 = disData;
    save('OD_DC_45_Degree_1.mat','data_OD_DC_45_Degree_1');
elseif strfind(fileName,'OD_DC_45_Degree_2')
    rotationMatrix = [cos(pi/4) -sin(pi/4) ; sin(pi/4) cos(pi/4)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    data_OD_DC_45_Degree_2 = disData;
    save('OD_DC_45_Degree_2.mat','data_OD_DC_45_Degree_2');
elseif strfind(fileName,'OD_DC_90_Degree_1') 
    rotationMatrix = [cos(-pi/2) -sin(-pi/2) ; sin(-pi/2) cos(-pi/2)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    data_OD_DC_90_Degree_1 = disData;
    save('OD_DC_90_Degree_1.mat','data_OD_DC_90_Degree_1');
elseif strfind(fileName,'OD_DC_90_Degree_2')
    disData(:,1) = -disData(:,1);
    data_OD_DC_90_Degree_2 = disData;
    save('OD_DC_90_Degree_2.mat','data_OD_DC_90_Degree_2');
elseif strfind(fileName,'OS_DC_0_Degree_1')
    disData(:,1) = -disData(:,1);
    data_OS_DC_0_Degree_1 = disData;
    save('OS_DC_0_Degree_1.mat','data_OS_DC_0_Degree_1');
elseif strfind(fileName,'OS_DC_0_Degree_2') 
    rotationMatrix = [cos(pi/2) -sin(pi/2) ; sin(pi/2) cos(pi/2)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    disData(:,1) = -disData(:,1);
    data_OS_DC_0_Degree_2 = disData;
    save('OS_DC_0_Degree_2.mat','data_OS_DC_0_Degree_2');
elseif strfind(fileName,'OS_DC_45_Degree_1')  
    rotationMatrix = [cos(-pi/4) -sin(-pi/4) ; sin(-pi/4) cos(-pi/4)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    disData(:,1) = -disData(:,1);
    data_OS_DC_45_Degree_1 = disData;
    save('OS_DC_45_Degree_1.mat','data_OS_DC_45_Degree_1');
elseif strfind(fileName,'OS_DC_45_Degree_2') 
    rotationMatrix = [cos(pi/4) -sin(pi/4) ; sin(pi/4) cos(pi/4)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    disData(:,1) = -disData(:,1);
    data_OS_DC_45_Degree_2 = disData;
    save('OS_DC_45_Degree_2.mat','data_OS_DC_45_Degree_2');
elseif strfind(fileName,'OS_DC_90_Degree_1')  
    rotationMatrix = [cos(-pi/2) -sin(-pi/2) ; sin(-pi/2) cos(-pi/2)];
    disData(:,1:2) = mtimes(rotationMatrix,disData(:,1:2)')';
    disData(:,1) = -disData(:,1);
    data_OS_DC_90_Degree_1 = disData;
    save('OS_DC_90_Degree_1.mat','data_OS_DC_90_Degree_1');
elseif strfind(fileName,'OS_DC_90_Degree_2')
    data_OS_DC_90_Degree_2 = disData;
    save('OS_DC_90_Degree_2.mat','data_OS_DC_90_Degree_2');
end

