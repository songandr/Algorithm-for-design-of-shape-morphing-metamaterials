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
% Updated Date: 07/15/25

% Initial x-values: planar kirigami, see Fig. S1 from "Programming bistability in
% geometrically perturbed mechanical metamaterials" Y. Peng et. al
%{
x1 = [0; 0];
x2 = [0.5; -1];
x3 = [0.4; 3];
x4 = [1.3; -1.1];
x5 = [2.6; 1.5];
x6 = [2; 4];
x7 = [2.2; 8];
x8 = [1.6; 7.5];
x9 = [0.7; 5];
x10 = x8+(x2-x4); % constrained for tessellation
x11 = x1+(x7-x5); % constrained for tessellation
x12 = [-0.3; 6];

l1R = x5-x1; % lattice vectors
l2R = x10-x2;
%}

% sheared planar kirigami, see Fig. 4 and Eq. (6) from above paper 
% gamma = 0.2; % Fig. 4 shows stability all the way up to gamma = 0.95
%l1D = l1R + gamma*l2R;
%l2D = l2R + gamma*l1R;
%{
x1 = [0; 0];
x2 = [1; 0];
x3 = [1; 1];
x4 = [1.5; 0.5];
x5 = [2; 1];
x6 = [1.8; 1.5];
x7 = [2.5; 2.5];
x8 = [2.2; 3];
x9 = [1.25; 1.25];
x10 = x8+(x2-x4); % constrained for tessellation
x11 = x1+(x7-x5); % constrained for tessellation
x12 = [0.3; 0.5];
%}

% Reference and deformed configuration lattice vector parameters
test = "shear";
%{
s_r = 0.5; % reference side length
l1R = 2*s_r*[1; 0]; % lattice vectors for fully rotated config
l2R = 2*s_r*[0; 1];

alpha1 = 0.5; % parameterization for cuts
alpha2 = 0.5;

x1 = [0; 0];
x5 = l1R;
x11 = l2R;
x7 = l1R+l2R;
x2 = x1 + alpha1*(x5-x1);
x4 = x2;
x8 = x11 + (1-alpha1)*(x7-x11);
x10 = x8;
x12 = x1 + alpha2*(x11-x1);
x6 = x5 + (1-alpha2)*(x7-x5);
x3 = (x2 + x8)/2;
x9 = x3;
%}

% Initial x-values
xi = pi/20; % pi/4 = fully rotated config
s = sqrt(2)/4; 
x1 = [0; 0];
x2 = s*[cos(xi); -sin(xi)];
x3 = s*[sin(xi)+cos(xi); cos(xi)-sin(xi)];
x4 = s*[2*sin(xi)+cos(xi); -sin(xi)];
x5 = s*[2*sin(xi)+2*cos(xi); 0];
x6 = s*[sin(xi)+2*cos(xi); cos(xi)];
x7 = s*[2*sin(xi)+2*cos(xi); 2*cos(xi)];
x8 = s*[2*sin(xi)+cos(xi); 2*cos(xi)+sin(xi)];
x9 = s*[sin(xi)+cos(xi); cos(xi)+sin(xi)];
x10 = s*[cos(xi); 2*cos(xi)+sin(xi)];
x11 = s*[0; 2*cos(xi)];
x12 = s*[sin(xi); cos(xi)];
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

l1R = x5 - x1;
l2R = x10 - x2;

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
% Initial x-values for planar kirigami initialized from lattice vectors & weighted averages
%{
x1 = [0; 0];
x5 = x1 + l1R;
x2 = 3/4 * x1 + 1/4 * x5;
x4 = 1/4 * x1 + 3/4 * x5;
x11 = x1 + l2R;
x12 = 1/2 * x11;
x7 = x5 + x11;
x6 = 1/2 * x5 + 1/2 * x7;
x10 = 3/4 * x11 + 1/4 * x7;
x8 = 1/4 * x11 + 3/4 * x7;
%x3 = 3/4 * 1/2 * (x6 + x12) + 1/4 * 1/2 * (x2 + x4);
%x9 = 3/4 * 1/2 * (x6 + x12) + 1/4 * 1/2 * (x8 + x10);
x3 = 3/4 * 1/2 * (x2 + x4) + 1/4 * 1/2 * (x10 + x8);
x9 = 3/4 * 1/2 * (x10 + x8) + 1/4 * 1/2 * (x2 + x4);
%}
%{
if test == "axial"
else
    y1 = [0; 0; 0];
    y2 = s*[cos(xi); -sin(xi); 0];
    y3 = s*[sin(xi)+cos(xi); cos(xi)-sin(xi); 0];
    y4 = s*[2*sin(xi)+cos(xi); -sin(xi); 0];
    y5 = s*[2*sin(xi)+2*cos(xi); 0; 0];
    y6 = s*[sin(xi)+2*cos(xi); cos(xi); 0];
    y7 = s*[2*sin(xi)+2*cos(xi); 2*cos(xi); 0];
    y8 = s*[2*sin(xi)+cos(xi); 2*cos(xi)+sin(xi); 0];
    y9 = s*[sin(xi)+cos(xi); cos(xi)+sin(xi); 0];
    y10 = s*[cos(xi); 2*cos(xi)+sin(xi); 0];
    y11 = s*[0; 2*cos(xi); 0];
    y12 = s*[sin(xi); cos(xi); 0];
end

%y = [y1; y2; y3; y4; y5; y6; y7; y8; y9; y10; y11; y12];
%}

% x and y rigidity constraint matrices (5x12) of (2x2) and (3x3)
L_0 = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2)];
L = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3)];
%{
U = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
A = [-eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];
%}

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% index set for panels 1~4 
F1 = [1, 2, 3, 12];
F2 = [3, 4, 5, 6];
F3 = [6, 7, 8, 9];
F4 = [9, 10, 11, 12];
T1 = F1;
T2 = F2;
T3 = F3;
T4 = F4;

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R
for j = 1:length(J)
    %R{j} = eye(3); % identity matrix
    if mod(j,2) == 1
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
tol = 10^(-5);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, L, L_0, tol);

% plotting
titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");

