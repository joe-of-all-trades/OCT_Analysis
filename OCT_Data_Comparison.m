% Thie script is written by Chao-Yuan Yeh. All rights reserved. 
%% Load and organize data into eight corners
clear
close all
if exist('OD_DC_0_Degree_1.mat','file'); load 'OD_DC_0_Degree_1.mat'; end
if exist('OD_DC_0_Degree_2.mat','file'); load 'OD_DC_0_Degree_2.mat'; end
load 'OD_DC_45_Degree_1.mat'
load 'OD_DC_45_Degree_2.mat'
if exist('OD_DC_90_Degree_1.mat','file')
    load 'OD_DC_90_Degree_1.mat'
    data_OD_DC_0_Degree_2 = data_OD_DC_90_Degree_1;
    clear data_OD_DC_90_Degree_1
end
if exist('OD_DC_90_Degree_2.mat','file')
    load 'OD_DC_90_Degree_2.mat'
    data_OD_DC_0_Degree_1 = data_OD_DC_90_Degree_2;
    data_OD_DC_0_Degree_1(:,2) = -data_OD_DC_0_Degree_1(:,2);
    clear data_OD_DC_90_Degree_2;
end
load 'OS_DC_0_Degree_1.mat'
load 'OS_DC_0_Degree_2.mat'
load 'OS_DC_45_Degree_1.mat'
load 'OS_DC_45_Degree_2.mat'

x1 = max(data_OD_DC_0_Degree_1(:,1));
x2 = abs(min(data_OD_DC_0_Degree_1(:,1)));
x3 = max(data_OD_DC_45_Degree_1(:,1))/sin(pi/4);
x4 = abs(min(data_OD_DC_45_Degree_1(:,1)))/sin(pi/4);
y1 = max(data_OD_DC_0_Degree_2(:,2));
y2 = abs(min(data_OD_DC_0_Degree_2(:,2)));
y3 = max(data_OD_DC_45_Degree_2(:,2))/sin(pi/4);
y4 = abs(min(data_OD_DC_45_Degree_2(:,2)))/sin(pi/4);
extent_OD = min([x1 x2 x3 x4 y1 y2 y3 y4]);

x1 = min(data_OD_DC_0_Degree_1(data_OD_DC_0_Degree_1(:,1)>0,1));
x2 = abs(max(data_OD_DC_0_Degree_1(data_OD_DC_0_Degree_1(:,1)<0,1)));
x3 = min(data_OD_DC_45_Degree_1(data_OD_DC_45_Degree_1(:,1)>0,1))/sin(pi/4);
x4 = abs(max(data_OD_DC_45_Degree_1(data_OD_DC_45_Degree_1(:,1)<0,1)))/sin(pi/4);
y1 = min(data_OD_DC_0_Degree_2(data_OD_DC_0_Degree_2(:,2)>0,2));
y2 = abs(max(data_OD_DC_0_Degree_2(data_OD_DC_0_Degree_2(:,2)<0,2)));
y3 = min(data_OD_DC_45_Degree_2(data_OD_DC_45_Degree_2(:,2)>0,2))/sin(pi/4);
y4 = abs(max(data_OD_DC_45_Degree_2(data_OD_DC_45_Degree_2(:,2)<0,2)))/sin(pi/4);
minextent_OD = max([x1 x2 x3 x4 y1 y2 y3 y4]);

x1 = max(data_OS_DC_0_Degree_1(:,1));
x2 = abs(min(data_OS_DC_0_Degree_1(:,1)));
x3 = max(data_OS_DC_45_Degree_1(:,1))/sin(pi/4);
x4 = abs(min(data_OS_DC_45_Degree_1(:,1)))/sin(pi/4);
y1 = max(data_OS_DC_0_Degree_2(:,2));
y2 = abs(min(data_OS_DC_0_Degree_2(:,2)));
y3 = max(data_OS_DC_45_Degree_2(:,2))/sin(pi/4);
y4 = abs(min(data_OS_DC_45_Degree_2(:,2)))/sin(pi/4);
extent_OS = min([x1 x2 x3 x4 y1 y2 y3 y4]);

x1 = min(data_OS_DC_0_Degree_1(data_OS_DC_0_Degree_1(:,1)>0,1));
x2 = abs(max(data_OS_DC_0_Degree_1(data_OS_DC_0_Degree_1(:,1)<0,1)));
x3 = min(data_OS_DC_45_Degree_1(data_OS_DC_45_Degree_1(:,1)>0,1))/sin(pi/4);
x4 = abs(max(data_OS_DC_45_Degree_1(data_OS_DC_45_Degree_1(:,1)<0,1)))/sin(pi/4);
y1 = min(data_OS_DC_0_Degree_2(data_OS_DC_0_Degree_2(:,2)>0,2));
y2 = abs(max(data_OS_DC_0_Degree_2(data_OS_DC_0_Degree_2(:,2)<0,2)));
y3 = min(data_OS_DC_45_Degree_2(data_OS_DC_45_Degree_2(:,2)>0,2))/sin(pi/4);
y4 = abs(max(data_OS_DC_45_Degree_2(data_OS_DC_45_Degree_2(:,2)<0,2)))/sin(pi/4);
minextent_OS = max([x1 x2 x3 x4 y1 y2 y3 y4]);

extent = min([extent_OD, extent_OS]);
minextent = max([minextent_OD minextent_OS]);

data_OD_DC_0_Degree_1 = data_OD_DC_0_Degree_1(abs(data_OD_DC_0_Degree_1(:,1)) >= minextent &...
    abs(data_OD_DC_0_Degree_1(:,1)) <= extent,:);
data_OD_DC_0_Degree_2 = data_OD_DC_0_Degree_2(abs(data_OD_DC_0_Degree_2(:,2)) >= minextent &...
    abs(data_OD_DC_0_Degree_2(:,2)) <= extent,:);
data_OD_DC_45_Degree_1 = data_OD_DC_45_Degree_1(abs(data_OD_DC_45_Degree_1(:,1)) >= minextent*sin(pi/4) &...
    abs(data_OD_DC_45_Degree_1(:,1))<=extent*sin(pi/4),:);
data_OD_DC_45_Degree_2 = data_OD_DC_45_Degree_2(abs(data_OD_DC_45_Degree_2(:,1)) >= minextent*sin(pi/4) &...
    abs(data_OD_DC_45_Degree_2(:,1))<=extent*sin(pi/4),:);

data_OS_DC_0_Degree_1 = data_OS_DC_0_Degree_1(abs(data_OS_DC_0_Degree_1(:,1)) >= minextent &...
    abs(data_OS_DC_0_Degree_1(:,1)) <= extent,:);
data_OS_DC_0_Degree_2 = data_OS_DC_0_Degree_2(abs(data_OS_DC_0_Degree_2(:,2)) >= minextent &...
    abs(data_OS_DC_0_Degree_2(:,2)) <= extent,:);
data_OS_DC_45_Degree_1 = data_OS_DC_45_Degree_1(abs(data_OS_DC_45_Degree_1(:,1)) >= minextent*sin(pi/4) &...
    abs(data_OS_DC_45_Degree_1(:,1))<=extent*sin(pi/4),:);
data_OS_DC_45_Degree_2 = data_OS_DC_45_Degree_2(abs(data_OS_DC_45_Degree_2(:,1)) >= minextent*sin(pi/4) &...
    abs(data_OS_DC_45_Degree_2(:,1))<=extent*sin(pi/4),:);

sumData_OD = cat(1, data_OD_DC_0_Degree_1, data_OD_DC_0_Degree_2,...
    data_OD_DC_45_Degree_1, data_OD_DC_45_Degree_2);

sumData_OS = cat(1,data_OS_DC_0_Degree_1,...
    data_OS_DC_0_Degree_2,data_OS_DC_45_Degree_1, data_OS_DC_45_Degree_2);
    
maxInd = max(max(sumData_OD(:,3)),max(sumData_OS(:,3)));

sumData_OD(:,3) = floor(sumData_OD(:,3)/maxInd*128);
colormap = jet(128);
sumData_OD(:,4:6) = squeeze(ind2rgb(sumData_OD(:,3),colormap));
figure;
scatter(sumData_OD(:,1),sumData_OD(:,2),[], sumData_OD(:,4:6)); axis equal

sumData_OS(:,3) = floor(sumData_OS(:,3)/maxInd*128);
colormap = jet(128);
sumData_OS(:,4:6) = squeeze(ind2rgb(sumData_OS(:,3),colormap));
figure;
scatter(sumData_OS(:,1),sumData_OS(:,2),[], sumData_OS(:,4:6)); axis equal
set(gca,'xdir','reverse')

%% save data
prompt ={'Enter patient number'};
dlg_title = 'Enter patient number';
num_lines = 1;
def = {''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
ptname = answer{1};
figure(1); 
title([ptname,' OD'])
saveas(gcf,[ptname,'_OD.jpg']);
figure(2); 
title([ptname,' OS'])
saveas(gcf,[ptname,'_OS.jpg']);

%% Calculate eight quadrants
OD_0_1 = datainterp(data_OD_DC_0_Degree_1,'x');
OD_0_1 = OD_0_1(2:end-1,:);
OD_0_2 = datainterp(data_OD_DC_0_Degree_2,'y');
OD_0_2 = OD_0_2(2:end-1,:);
OD_45_1 = datainterp(data_OD_DC_45_Degree_1,'45');
OD_45_2 = datainterp(data_OD_DC_45_Degree_2,'45');

OD_Q1 = OD_0_1(OD_0_1(:,1)>0,:);
OD_Q2 = OD_45_2(OD_45_2(:,1)>0 & OD_45_2(:,2)>0, :);
OD_Q3 = OD_0_2(OD_0_2(:,2)>0,:);
OD_Q4 = OD_45_1(OD_45_1(:,1)<0 & OD_45_1(:,2)>0, :);
OD_Q5 = OD_0_1(OD_0_1(:,1)<0,:);
OD_Q6 = OD_45_2(OD_45_2(:,1)<0 & OD_45_2(:,2)<0, :);
OD_Q7 = OD_0_2(OD_0_2(:,2)<0,:);
OD_Q8 = OD_45_1(OD_45_1(:,1)>0 & OD_45_1(:,2)<0, :);

OS_0_1 = datainterp(data_OS_DC_0_Degree_1,'x');
OS_0_1 = OS_0_1(2:end-1,:);
OS_0_2 = datainterp(data_OS_DC_0_Degree_2,'y');
OS_0_2 = OS_0_2(2:end-1,:);
OS_45_1 = datainterp(data_OS_DC_45_Degree_1,'45');
OS_45_2 = datainterp(data_OS_DC_45_Degree_2,'45');

OS_Q1 = OS_0_1(OS_0_1(:,1)>0,:);
OS_Q2 = OS_45_1(OS_45_1(:,1)>0 & OS_45_1(:,2)>0, :);
OS_Q3 = OS_0_2(OS_0_2(:,2)>0,:);
OS_Q4 = OS_45_2(OS_45_2(:,1)<0 & OS_45_2(:,2)>0, :);
OS_Q5 = OS_0_1(OS_0_1(:,1)<0,:);
OS_Q6 = OS_45_1(OS_45_1(:,1)<0 & OS_45_1(:,2)<0, :);
OS_Q7 = OS_0_2(OS_0_2(:,2)<0,:);
OS_Q8 = OS_45_2(OS_45_2(:,1)>0 & OS_45_2(:,2)<0, :);

%% Save Statistics
xlsdata = {'OD_Q1',mean(OD_Q1(:,3));'OD_Q2', mean(OD_Q2(:,3));'OD_Q3', mean(OD_Q3(:,3));...
    'OD_Q4', mean(OD_Q4(:,3)); 'OD_Q5', mean(OD_Q5(:,3)); 'OD_Q6', mean(OD_Q6(:,3));...
    'OD_Q7', mean(OD_Q7(:,3)); 'OD_Q8', mean(OD_Q8(:,3));...
    'OD_Avg', mean(cat(1,OD_Q1(:,3),OD_Q2(:,3),OD_Q3(:,3),OD_Q4(:,3),OD_Q5(:,3),OD_Q6(:,3),OD_Q7(:,3),OD_Q8(:,3)));...
    'OD_Avg_Q4-6', mean([mean(OD_Q4(:,3)) mean(OD_Q5(:,3)) mean(OD_Q6(:,3))]);'',[];
    'OS_Q1',mean(OS_Q1(:,3)); 'OS_Q2', mean(OS_Q2(:,3));...
    'OS_Q3', mean(OS_Q3(:,3));'OS_Q4', mean(OS_Q4(:,3));...
    'OS_Q5', mean(OS_Q5(:,3));'OS_Q6', mean(OS_Q6(:,3));...
    'OS_Q7', mean(OS_Q7(:,3));'OS_Q8', mean(OS_Q8(:,3));...
    'OS_Avg', mean(cat(1,OS_Q1(:,3),OS_Q2(:,3),OS_Q3(:,3),OS_Q4(:,3),OS_Q5(:,3),OS_Q6(:,3),OS_Q7(:,3),OS_Q8(:,3)));...
    'OS_Avg_Q4-6', mean([mean(OS_Q4(:,3)) mean(OS_Q5(:,3)) mean(OS_Q6(:,3))])};
xlswrite([ptname,'_stat.xls'],xlsdata)

%% Smooth data
OD_Q1 = clean_and_smooth_data(OD_Q1);
OD_Q2 = clean_and_smooth_data(OD_Q2);
OD_Q3 = clean_and_smooth_data(OD_Q3);
OD_Q4 = clean_and_smooth_data(OD_Q4);
OD_Q5 = clean_and_smooth_data(OD_Q5);
OD_Q6 = clean_and_smooth_data(OD_Q6);
OD_Q7 = clean_and_smooth_data(OD_Q7);
OD_Q8 = clean_and_smooth_data(OD_Q8);

OS_Q1 = clean_and_smooth_data(OS_Q1);
OS_Q2 = clean_and_smooth_data(OS_Q2);
OS_Q3 = clean_and_smooth_data(OS_Q3);
OS_Q4 = clean_and_smooth_data(OS_Q4);
OS_Q5 = clean_and_smooth_data(OS_Q5);
OS_Q6 = clean_and_smooth_data(OS_Q6);
OS_Q7 = clean_and_smooth_data(OS_Q7);
OS_Q8 = clean_and_smooth_data(OS_Q8);

%% Interpret into circle

colormap = jet(1000);
template = generate_circle({OD_Q1, OD_Q2, OD_Q3, OD_Q4, OD_Q5, OD_Q6, ...
  OD_Q7,OD_Q8}, maxInd);
image_OD = squeeze(ind2rgb(template,colormap));
figure(3); imshow(image_OD)

template = generate_circle({OS_Q1, OS_Q2, OS_Q3, OS_Q4, OS_Q5, OS_Q6, ...
  OS_Q7, OS_Q8}, maxInd);
image_OS = squeeze(ind2rgb(template,colormap));
image_OS = flip(image_OS,2);
figure(4); imshow(image_OS)

imwrite(image_OD, [ptname,' OD_full.tif']);
imwrite(image_OS, [ptname,' OS_full.tif']);


