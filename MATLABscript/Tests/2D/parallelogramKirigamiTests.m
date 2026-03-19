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
% Updated Date: 01/26/26.

% Initial x-values: planar kirigami, see Fig. 2 from "Mechanism-based
% metamaterials with microstructurally invariant shape-change" Y. Peng et. al

% Initial x-values
x = []; % see unit_cell_design_recipe.m

l1R = x21 - x1;
l2R = x16 + (x21-x20) - x1;

phi = zeros(length(x)/2*3, length(x));
if test == "reference"
    %{
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
        end
    end
    %}
    y1 = [0; 0; 0];
    y5 = [l1R(1); l1R(2); 0] + y1;
    y2 = (y1+y5)/2;
    y4 = y2;
    y10 = [l2R(1); l2R(2); 0] + y2;
    y8 = y10;
    y11 = [l2R(1); l2R(2); 0] + y1;
    y12 = (y1+y11)/2;
    y7 = [l2R(1)+l1R(1); l2R(2)+l1R(2); 0] + y1;
    y6 = (y5+y7)/2;
    y3 = (y1+y7)/2;
    y9 = y3;
    y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12];
elseif test == "axial"
    lambda_1 = 1.3; % axial deformation
    lambda_2 = 1.3;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end
    y = phi*x;
elseif test == "shear"
    x1_box = x1;
    x2_box = l1R/2;
    x4_box = x2_box;
    x5_box = l1R;
    x6_box = x5_box+l2R/2;
    x7_box = x5_box+l2R;
    x11_box = l2R;
    x12_box = l2R/2;
    x3_box = x12_box + l1R/2;
    x9_box = x3_box;
    x8_box = x11_box + l1R/2;
    x10_box = x8_box;
    x_box = [x1_box; x2_box; x3_box; x4_box; x5_box; x6_box; x7_box; x8_box; x9_box; x10_box; x11_box; x12_box];
    
    gamma = 0.90;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)+1) = gamma;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
            phi(i,i-floor(i/3)-1) = gamma;
        end    
    end
    y = phi*x_box;
end

% x and y rigidity constraint matrices (12x23) of (2x2) and (3x3)
L_0 = zeros(2*12,2*23);
L_0(1:2, 1:2) = eye(2); % fix x1
L_0(3:4, 1:2) = -eye(2);
L_0(3:4, 21*2-1:21*2) = eye(2); % x21 - x1 = l1R
L_0(5:6, 4*2-1:4*2) = -eye(2);
L_0(5:6, 22*2-1:22*2) = eye(2); % x22 - x4 = l1R
L_0(7:8, 9*2-1:9*2) = -eye(2);
L_0(7:8, 23*2-1:23*2) = eye(2); % x23 - x9 = l1R
L_0(9:10, 11*2-1:11*2) = -eye(2);
L_0(9:10, 14*2-1:14*2) = eye(2); 
L_0(9:10, 23*2-1:23*2) = eye(2);
L_0(9:10, 13*2-1:13*2) = -eye(2); % x14 + (x23-x13) - x11 = l1R
L_0(11:12, 17*2-1:17*2) = -eye(2);
L_0(11:12, 21*2-1:21*2) = eye(2);
L_0(11:12, )
L_0(13:14, 1:2) = -eye(2);
L_0(13:14, 16*2-1:16*2) = eye(2);
L_0(13:14, 21*2-1:21*2) = eye(2);
L_0(13:14, 20*2-1:20*2) = -eye(2); % x16 + (x21-x20) - x1 = l2R
L_0(15:16, 17*2-1:17*2) = -eye(2);


L = zeros(3*9,3*23);
% populate the vector numbering all of the panels
J = [1, 2, 3, 4, 5, 6, 7, 8];

% index set for panels 1~4 
F1 = [1, 2, 3, 4];
F2 = [5, 6, 7, 8];
F3 = [9, 10, 11];
F4 = [12, 13, 14, 15, 16];
F5 = [2, 5, 17];
F6 = [6, 18, 19, 20, 21];
F7 = [3, 8, 10, 17];
F8 = [7, 13, 22, 23];
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;
T5 = F5;
T6 = F6;
T7 = F7;
T8 = F8;

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R
for j = 1:length(J)
    %R{j} = eye(3); % identity matrix
    if j < 5
        epsilon = 0.0001;
    else
        epsilon = -0.0001;
    end
    W = [0  -1  0;
         1  0   0;
         0  0   0];
    R{j} = eye(3) + epsilon * W;
end

% initial tolerance for minimization
tol = 10^(-4);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, L_0, tol);

% plotting
titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");

