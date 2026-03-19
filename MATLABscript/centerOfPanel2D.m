function [cj, rij] = centerOfPanel2D(Tc, coord)
% Inputs:
% Tc: The set of all of the indices (as opposed to position) of x's within a panel
% coord: 2D array of all of the x coordinates (2*n by 1 where n is the number of indices)
%
% Outputs:
% cj: the center of the panel
% rij: 2D array containing the coordinate vectors with respect to the center of the panel
% ordered sequentially so iterating through i works 

len = length(Tc);
sum = zeros(2, 1);

for i = 1:len
    sum = sum + coord(2*Tc(i)-1:2*Tc(i), 1);
end

cj = sum/len;

rij = zeros(2*len, 1);
for i = 1:len
    rij(2*i-1:2*i, 1) = coord(2*Tc(i)-1:2*Tc(i), 1) - cj;
end


end

