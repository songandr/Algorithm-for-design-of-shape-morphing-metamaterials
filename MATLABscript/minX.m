function xMin = minX(x_c, y, Fj, J, R, L_0, A_0, xMap)
%
% directly solves the linear equation to calculate the perturbation of the
% initial x
%
% Inputs:
%
% Rigidity Constraints:
% L_0: matrix of the rigidity constraints that satisfies the equation: L_0*x = c  
%
% Indexing Inputs:
% x_c: x coordinate 2-D array (2*n by 1 where n is the number of vertices)
% y: y coordinate 2-D array (3*n by 1 where n is the number of vertices)
% Fj: a cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% J: the labeling set of all panels
% R: a cell array of all of the rotation matrices for each panel
% A_0: matrix pre-computed for vertex minimization steps
% xMap: cell array storing chihat_i - avg(chihat)_j in {i, j}th entry
%
% Outputs
% xMin: x coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   

% Initialize relevant parameters
n = length(x_c); % number of vertices * 2
lenJ = length(J); % number of panels

% determine null space of L_0
N_T = null(L_0);

% pos vectors with respect to the center of the panel for all panels
rij = cell(lenJ); % j cells with the j-th cell containing 
for j = 1:lenJ
    [~, rij{j}] = centerOfPanel3D(Fj{j}, y);
end

a_0 = zeros(n, 1);
for j = 1:lenJ
    for i = 1:length(Fj{j})
       k = Fj{j}(i); % vertex k
       a_0 = a_0 + xMap{k,j}'*R{j}'*rij{j}(3*i-2:3*i);
    end
end

% calculation of aTilde
aTilde = N_T'*a_0 - N_T'*A_0*x_c ; 

% determining the perturbation method
aTilde(abs(aTilde)<1e-5)=0;

z_0 = (N_T'*A_0*N_T)\aTilde;

xMin = x_c + N_T*z_0;

end