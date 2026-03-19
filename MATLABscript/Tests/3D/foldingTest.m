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
% Updated Date: 12/06/24
%
% The initial configuration consists of two rectangular panels. This test
% is interested in the folding of these two panels in the z-axis.

% Initial x-values
x1 = [0; 0];
x2 = [1; 0];
x3 = [2; 0];
x4 = [2; 1];
x5 = [1; 1];
x6 = [0; 1];
x = [x1; x2; x3; x4; x5; x6];

% Initial y-values: ICs such that there is folding upwards in z's for
% resultant lattice vectors in the form (lambda, 0) where lambda < 1
lambda = 0.5;
y1 = [0; 0; 0];
y2 = [lambda; 0; sqrt(1-lambda^2)+3];
y3 = [2*lambda; 0; 0];
y4 = [2*lambda; 1; 0];
y5 = [lambda; 1; sqrt(1-lambda^2)+3];
y6 = [0; 1; 0];
y = [y1; y2; y3; y4; y5; y6];

% x rigidity constraint matrix (6x6) of (2x2) = (12x12)
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);
     zeros(2), -eye(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), eye(2), zeros(2), zeros(2);];

% y rigidity constraint matrix (6x6) of (3x3) = (18x18)
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);
     zeros(3), -eye(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), eye(3), zeros(3), zeros(3);];

J = [1, 2];
F1 = [1, 2, 5, 6];
F2 = [2, 3, 4, 5];
T1 = F1;
T2 = F2;
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Fj = Tj;

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for minimization
tol = 10^(-3);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Folding");