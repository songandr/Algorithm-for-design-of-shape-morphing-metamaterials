function yMin = minY(x, y_d, Tj, J, R, L, A, yMap)
%
% directly solves the linear equation to calculate the perturbation of the
% initial y 
%
% Inputs
%
% Rigidity Constraints:
% L: matrix of the rigidity constraints that satisfies the equation: L*y = d  
%
% Indexing Inputs:
% x: x coordinate 2-D array (2*n by 1 where n is the number of vertices)
% y_d: y coordinate 2-D array (3*n by 1 where n is the number of vertices)
% Tj: a cell array of the set of all x's within each panel (the jth panel
% corresponds to the jth row)
% J: the labeling set of all panels
% R: a cell array of all of the rotation matrices for each panel
% A: matrix pre-computed for vertex minimization steps
% yMap: cell array storing chi_i - avg(chi)_j in {i, j}th entry
%
% Outputs
% yMin: y coordinate 2-D array that minizes the elastic energy based on
% given rigidity constraints   

% Initialize relevant parameters
n = length(y_d); % number of vertices * 3
lenJ = length(J); % number of panels

% determine null space of L
N_G = null(L);

% pos vectors with respect to the center of the panel for all panels
rij = cell(lenJ); % j cells with the j-th cell containing 
for j = 1:lenJ
    [~, rij{j}] = centerOfPanel2D(Tj{j}, x);
end 

a = zeros(n, 1);
for j = 1:lenJ
    for i = 1:length(Tj{j})
        k = Tj{j}(i); % vertex k
        r_temp = rij{j}(2*i-1:2*i);
        rij1 = r_temp(1);
        rij2 = r_temp(2);
        a = a + yMap{k,j}'*R{j}*[rij1; rij2; 0]; % see Eq. 17
    end
end

% calculation of bTilde
bTilde = N_G'*a - N_G'*A*y_d;

% determining the perturbation method
bTilde(abs(bTilde)<1e-5)=0;

% condition matrix for inverting
preinverse_tensor = N_G'*A*N_G;
preinverse_tensor = (preinverse_tensor+preinverse_tensor')/2; % symmetricize tensor before inverting
preinverse_tensor(abs(preinverse_tensor)<1e-10) = 0; % round tensor before inverting
z = preinverse_tensor\bTilde;
disp("rcond: "+rcond(N_G'*A*N_G))

yMin = y_d + N_G*z;

end