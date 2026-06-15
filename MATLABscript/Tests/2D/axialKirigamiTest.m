function [xOpt, yOpt, E] = axialKirigamiTest(x, y, L_0, L)
% 
% runs 1 axial shape change test on a rotating squares topology with N vertices
%
% inputs:
% x, y: a 2N and 3N sized column array holding reference and deformed configuration vertices
% L, L_0: x and y rigidity constraint matrices (5x12) of (2x2) and (3x3)
%
% outputs:
% xOpt: energy minimizing reference configuration vertices (length 2N)
% yOpt: energy minimizing deformed configuration vertices (length 3N)

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% index set for panels 1~4 
P1 = [1, 2, 3, 12];
P2 = [3, 4, 5, 6];
P3 = [6, 7, 8, 9];
P4 = [9, 10, 11, 12];

% 3D array containing index set of x coordinates for panel j
Pj = cell(length(J), 1);
for j=1:length(J)
    Pj{j} = eval(sprintf('P%d', j));
end

% Initial R 
R = cell(1, length(J));
for j = 1:length(J)
    R{j} = eye(3); % identity
end

% tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Pj, J, R, L, L_0, tol);

% Compute energy associated with solution
E = 0;

% Construct necessary cj and rij vectors
cj = cell(length(J));
rij = cell(length(J));
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel3D(Pj{j}, yOpt);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Pj{j}, xOpt);
end

for j = 1:length(J)
    for i = 1:length(Pj{j})
        k = Pj{j}(i);

        rij_temp = rij{j}(2*i-1:2*i);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);
        
        E = E + norm(yOpt(3*k-2:3*k, 1) - cj{j} - Ropt{j}*[rij1; rij2; 0])^2;
    end
end 
end