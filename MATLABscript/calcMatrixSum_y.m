function matrixSum = calcMatrixSum_y(chi, Fj)
% 
% calculates the summation of matrices chi that selects the y_i vector
%
% inputs:
% Fj: an array of the set of all y within the given panel 
% xM: cell array of matrices
%
% outputs:
% matrixSum: the sum of all of the matrices corresponding to a given panel
%

n = length(chi);
sum = zeros(3, 3*n);

for k = 1:length(Fj)
   p = Fj(k);
   sum = sum + chi{p};
end

matrixSum = sum;

end