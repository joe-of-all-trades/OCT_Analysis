function output = clean_and_smooth_data( input )
% Use average filter to smooth data. 
kernel = [.2 .2 .2 .2 .2];
input(:,3) = filter(kernel, 1, input(:,3));
input(end-2:end,:) = [];
output = check_order(input);
end

