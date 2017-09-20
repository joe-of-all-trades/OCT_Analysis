function output = generate_circle( input, maxInd )
% Use data in eight corners and interpolate into a circle.
% This script is written by Chao-Yuan Yeh.
% All rights reserved. 

output = zeros(800,800);

for ii = 1:min([length(input{1}) length(input{2}) length(input{3}) ...
        length(input{4}) length(input{5}) length(input{6}) length(input{8})])
    x = pi*(0:.25:2);
    y = input{1}(ii,1)*[0  1 cos(pi/4) 0 -cos(pi/4) -1 -cos(pi/4) 0 cos(pi/4) 1 0;
        1  0 sin(pi/4) 1 sin(pi/4) 0 -sin(pi/4) -1  -sin(pi/4) 0 1];
    pp = spline(x,y);
    bq = 2*pi*input{1}(ii,1);
    yy = ppval(pp, linspace(0,2*pi,bq))';
    yy = yy(2:end-1,:);
    
    basedata = [input{1}(ii,3) input{2}(ii,3) input{3}(ii,3) input{4}(ii,3)...
        input{5}(ii,3) input{6}(ii,3) input{7}(ii,3) input{8}(ii,3)...
        input{1}(ii,3) input{2}(ii,3) input{3}(ii,3)];
    
    base = [0 floor(bq/8) floor(bq/4) floor(bq*3/8) floor(bq/2)...
        floor(bq*5/8) floor(bq*3/4) floor(bq*7/8) bq floor(9*bq/8) floor(5*bq/4)];

    bqs = 1:floor(5*bq/4);
    yy2 = interp1(base,basedata,bqs,'spline')';
    yy2 = yy2(floor(bq/8:floor(bq/8)+length(yy)));
    yy2 = circshift(yy2,floor(bq/8));
    yy(:,3) = round(yy2);
    
    xshift = abs(floor(min(yy(:,1))));
    yshift = abs(floor(min(yy(:,2))));
    
    yy(:,1) = floor(yy(:,1))+xshift+1+floor((800-range(yy(:,1)))/2);
    yy(:,2) = floor(yy(:,2))+yshift+1+floor((800-range(yy(:,2)))/2);
    
    gap = max(abs(yy(1,1)-yy(end,1)), abs(yy(end,2)-yy(1,2)))-1;
    gapfill = round((yy(end,3) + yy(1,3))/2);
    for jj = 1:gap
        yy(end+1,:) = [ yy(end,1)+(yy(end,1)-yy(1,1)~=0),...
            yy(end,2) + (yy(end,2)-yy(1,2)~=0), gapfill];
    end
    
    for jj = 1:length(yy)
        output(yy(jj,2),yy(jj,1)) = yy(jj,3);
    end
end

output = medfilt2(output);
kernel = fspecial('gaussian',7,.7);
output = imfilter(output,kernel);
output = floor(output/maxInd*1000);
output = flip(output,1);
end

