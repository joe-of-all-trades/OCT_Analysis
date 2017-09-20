function output = check_order( input )
% Make sure the list is in descending order in terms of distance from
% origin.
if sum(input(1, 1:2).^2) > sum(input(end, 1:2).^2)
  output = input;
else
  output = flip(input, 1);
end

end