function [chi, chihat] = construct_chi(y)
% 
% constructs the maps chi and chihat s.t. chi_i*y = y_i, chihat_i*x = (x_i, 0)
%
% Inputs
% y: y coordinate 2-D column array (3*n by 1 where n is the number of vertices)
% 
% Outputs
% chi: cell array where i-th entry is the 3 by 3*n mapping matrix mapping y
% in R^(3n) to y_i in R^3
% chihat: cell array where i-th entry is the 3 by 2*n mapping matrix
% mapping x in R^(2n) to (x_i, 0) in R^3

% Initialize relevant parameters
n = length(y); % number of vertices * 3
chi = cell(n/3,1);
chihat = cell(n/3,1);

for k = 1:(n/3)
    chi_k = zeros(3, n);
    chi_k(1:3, 3*k-2:3*k) = eye(3);
    chi{k} = chi_k;

    chihat_k = zeros(3, 2*n/3);
    chihat_k(1:3, 2*k-1:2*k) = [1  0;
                                0  1;
                                0  0];
    chihat{k} = chihat_k;
end