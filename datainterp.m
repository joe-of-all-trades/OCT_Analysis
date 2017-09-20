% This script is written by Chao-yuan Yeh. All copyrights reserved. 
function output = datainterp( input, varargin)

if strcmpi(varargin, '45')
    temp = diff(round(input(:,1)/sin(pi/4)));
elseif strcmpi(varargin, 'x')
    temp = diff(round(input(:,1)));
elseif strcmpi(varargin, 'y')
    temp = diff(round(input(:,2)));
end

index = find(abs(temp)==2);
index = index + (1:length(index))';
temp2 = zeros(size(input,1)+length(index), size(input,2));
temp2(index,:) = nan;
temp2(~isnan(temp2)) = input;
temp2(index,:) = (temp2(index-1,:) + temp2(index+1,:))/2;
output = temp2;

end

