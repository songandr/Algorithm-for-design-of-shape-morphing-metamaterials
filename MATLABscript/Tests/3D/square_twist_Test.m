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
% Updated Date: 04/09/25.
%
% The initial configuration consists of nine panels in the square twist
% configuration. This test is interested in the folding of these nine 
% panels and the resultant energy calculation.

% Initial x-values
s = 1;
phi = pi/6; % rhombus acute angle
x1 = [0; 0];
x2 = s*[1; 0];
x3 = s*[1+sin(phi); -cos(phi)];
x4 = s*[2+sin(phi); -cos(phi)];
x5 = s*[2+sin(phi); -cos(phi)+1];
x6 = s*[1+sin(phi); -cos(phi)+1];
x7 = s*[1; 1];
x8 = s*[0; 1];
x9 = s*[cos(phi); 1+sin(phi)];
x10 = s*[cos(phi)+1; 1+sin(phi)];
x11 = s*[1+sin(phi)+cos(phi); -cos(phi)+1+sin(phi)];
x12 = s*[2+sin(phi)+cos(phi); -cos(phi)+1+sin(phi)];
x13 = s*[2+sin(phi)+cos(phi); -cos(phi)+2+sin(phi)];
x14 = s*[1+sin(phi)+cos(phi); -cos(phi)+2+sin(phi)];
x15 = s*[cos(phi)+1; 2+sin(phi)];
x16 = s*[cos(phi); 2+sin(phi)];
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13; x14; x15; x16];
l1R = x4 - x1;
l2R = x16 - x1;

% Initial y-values
theta = phi + pi/4; % twist angle
alpha = pi/180; % how low is the center square relative to outer squares
y1 = s*[0; 0; 0];
y2 = s*[1; 0; 0];
y3 = s*[1-cos(pi/2 - theta + phi)*cos(alpha); -sin(pi/2 - theta + phi)*cos(alpha); -sin(alpha)]; % from y4
y4 = s*[2-cos(pi/2 - theta + phi)*cos(alpha); -sin(pi/2 - theta + phi)*cos(alpha); -sin(alpha)]; % from y5
y5 = s*[2-cos(pi/2 - theta + phi)*cos(alpha); 1-sin(pi/2 - theta + phi)*cos(alpha); -sin(alpha)]; % from y6
y6 = s*[1-cos(pi/2 - theta + phi)*cos(alpha); 1-sin(pi/2 - theta + phi)*cos(alpha); -sin(alpha)]; % from y7
y7 = s*[1; 1; 0];
y8 = s*[0; 1; 0];
y10 = y7 + [0 -1 0; 1 0 0; 0 0 1] * (y6-y7); % from y7
y9 = y10 - s*[1; 0; 0]; % from y10
y11 = y10 + [0 -1 0; 1 0 0; 0 0 1] * (y7-y10); % from y10
y12 = y11 + s*[1; 0; 0]; % from y11
y13 = y12 + s*[0; 1; 0]; % from y12
y14 = y13 - s*[1; 0; 0]; % from y13
y15 = y10 + s*[0; 1; 0]; % from y10
y16 = y15 - s*[1; 0; 0]; % from y15
y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12; y13; y14; y15; y16];

% x rigidity constraint matrix (4x16) of (2x2) = (8x32)
U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2)];
    
% x rigidity constraint matrix (4x12) of (3x3) = (12x36)
A = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), -eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3)];
% populate the vector numbering all of the panels
J = [1, 2, 3, 4, 5, 6, 7, 8, 9];

% determine the index set for each panel 
% a cell array of the set of all y's within each panel
F1 = [1, 2, 7, 8];
F2 = [2, 3, 6, 7];
F3 = [3, 4, 5, 6];
F4 = [5, 6, 11, 12];
F5 = [6, 7, 10, 11];
F6 = [7, 8, 9, 10];
F7 = [9, 10, 15, 16];
F8 = [10, 11, 14, 15];
F9 = [11, 12, 13, 14];

% a cell array of the set of all x's within each panel 
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;
T5 = F5;
T6 = F6;
T7 = F7;
T8 = F8;
T9 = F9;

% a 3-d array that contains the index set of x coordinates within a given
% panel j
Tj = zeros(1, length(T1), length(J));
Tj(1, :, 1) = T1;
Tj(1, :, 2) = T2;
Tj(1, :, 3) = T3;
Tj(1, :, 4) = T4;
Tj(1, :, 5) = T5;
Tj(1, :, 6) = T6;
Tj(1, :, 7) = T7;
Tj(1, :, 8) = T8;
Tj(1, :, 9) = T9;

Fj = Tj;

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% determine the initial tolerance for minimization
tol = 10^(-4);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, "Square Twist");
