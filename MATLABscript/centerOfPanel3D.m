function [cj, rij] = centerOfPanel3D(Fc, coord)
% Inputs:
% Fc: The set of all of the indices (as opposed to position) of y's within a panel
% coord: 2D array of all of the y coordinates (3*n by 1 where n is the number of indices)
%
% Outputs:
% cj: the center of the panel
% rij: 2D array containing the coordinate vectors with respect to the center of the panel
% ordered sequentially so iterating through i works 

len = length(Fc);
sum = zeros(3, 1);

for i = 1:len
    sum = sum + coord(3*Fc(i)-2:3*Fc(i), 1);
end

cj = sum/len;

rij = zeros(3*len, 1);
for i = 1:len
    rij(3*i-2:3*i, 1) = coord(3*Fc(i)-2:3*Fc(i), 1) - cj;
end


end

