function matrixSum = calcMatrixSum_x(chihat, Pj)
% 
% calculates the summation of matrices chihat that selects the x_i vector
%
% inputs:
% Pj: an array of the set of all x within a given panel 
% chihat: cell array of mapping matrices
%
% outputs:
% matrixSum: the sum of all of the matrices corresponding to a given panel

n = length(chihat);
sum = zeros(3, 2*n);

for k = 1:length(Pj)
   p = Pj(k);
   sum = sum + chihat{p};
end

matrixSum = sum; 

end

