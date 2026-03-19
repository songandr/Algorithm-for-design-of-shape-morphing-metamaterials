% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Antoine Moats, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 07/29/25.
%
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.

% Initial x-values
s = 0.5;
theta_0 = pi/2*0.95;
mu_0 = 1;

l1R = [2*s; 0];
l2R = mu_0 * [cos(theta_0)  -sin(theta_0); sin(theta_0) cos(theta_0)]*l1R;

x1 = [0; 0];
x3 = l1R;
x7 = l2R;
x6 = x7/2;
x4 = x6 + l1R;
x9 = x7 + l1R;
x2 = (x1+x3+x4+x6)/4;
x8 = x2 + l2R;
x5 = (x2+x8)/2;
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

lambda = 0.9;
theta = pi/4;
mu = mu_0 * cos(theta_0) / (lambda * cos(theta)); % enforce rigidity

l1D = [lambda; 0; 0];
l2D = mu*[cos(theta)  -sin(theta)   0; sin(theta) cos(theta)    0;  0   0   1]*[1; 0; 0];

y1 = [0; 0; 0];
y3 = l1D;
y7 = l2D;
y9 = l1D + l2D;
y4 = (y3+y9)/2;
y6 = y7/2;
y2 = (y1+y3+y4+y6)/4;
y8 = y2 + l2D;
y5 = (y2+y8)/2;
% Bias y's (y4 y5 y6) in the z-axis
y4(3) = y4(3)+0.1;
y5(3) = y5(3)+0.1;
y6(3) = y6(3)+0.1;

y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];

% x rigidity constraint matrix (7x9) of (2x2)
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);];
    
% y rigidity constraint matrix (7x9) of (3x3)
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);];

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% index set for each panel 
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
F3 = [5, 6, 7, 8];
F4 = [4, 5, 8, 9];
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;

% 3D array containing index set of x coordinates for panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Tj(1, :, 3) = T3;
Tj(1, :, 4) = T4;
Fj = Tj;

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");
