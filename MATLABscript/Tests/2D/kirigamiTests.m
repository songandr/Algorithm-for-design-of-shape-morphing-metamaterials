% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Algorithmic design framework for shape-morphing metamaterials"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 06/15/26.
close all;
clear all; 
clc;

% Reference and deformed configuration lattice vector parameters
test = "axial";
star = 0; % for star-shaped slits
trapezoid = 1; % for trapezoid-shaped slits
check_strain = 0;

% Initial x-values
%xi = pi/4; % for fully open config
xi = 0.1; % for negligibly opened config
%s = sqrt(2)/4; % for fully open config
s = 0.5*0.9134; % for negligibly opened config
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

l1R = x5 - x1;
l2R = x10 - x2;

if star
l1R = [1; 0];
l2R = [0; 1];
x1 = [0; 0];
x11 = x1 + l2R;
x12 = (x1+x11)/2;
x6 = x12 + l1R;
x5 = x1 + l1R;
x7 = x5 + l2R;
x2 = x1 + l1R/4;
x4 = x5 - l1R/4;
x3 = x1 + l1R/2 + l2R/4;
x10 = x11 + l1R/4;
x8 = x7 - l1R/4;
x9 = x11 + l1R/2 - l2R/4;
elseif trapezoid
s = 0.5; % side length
b = -0.25; % departure from parallelogram
l1R = s*[2; 2];
l2R = s*[-2; 2];
x1 = [0; 0]; % P1 square
x2 = [s; 0];
x3 = [s; s];
x12 = [0; s];
x6 = x3 + [b; s]; % construct trapezoid
x9 = x12 + [0; s];
x4 = [2*s; s]; % P2
x5 = x1 + l1R;
x10 = x2 + l2R; % P4
x11 = [-s; s]; 
x7 = x11 + l1R; % P3
x8 = x4 + l2R;
end

x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

% perturb initialization to avoid getting stuck
if trapezoid
    lambda = norm(x6-x9)/norm(x3-x12); % trapezoid parallel side scaling
    rng(1, "twister");
    r = normrnd(0, 0.05, [12, 1]); % only 6 DoFs to maintain original lattice vectors post-perturbation
    x(3:4) = x(3:4)+r(1:2); % 1. x2
    x(19:20) = x(19:20)+r(1:2); % x10 follows x2
    x(5:6) = x(5:6)+r(3:4); % 2. x3
    x(7:8) = x(7:8)+r(5:6); % 3. x4
    x(15:16) = x(15:16)+r(5:6); % x8 follows x4
    x(13:14) = x(13:14)+r(7:8); % 4. x7
    x(21:22) = x(21:22)+r(7:8); % x11 follows x7
    x(17:18) = x(17:18)+r(9:10); % 5. x9
    x(23:24) = x(23:24)+r(11:12); % 6. x12
    x(11:12) = x(11:12)+lambda*(r(3:4)-r(11:12))+r(9:10); % x6 must satisfy trapezoid constraint
end

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
    lambda_1 = 0.7; % axial deformation
    lambda_2 = 1;
    if trapezoid % spectral decomposition along diagonals
        for i = 1:length(x)/2*3
            if mod(i,3) == 1
                phi(i,i-floor(i/3)) = lambda_1+lambda_2;
                phi(i,i-floor(i/3)+1) = lambda_1-lambda_2;
            elseif mod(i,3) == 2
                phi(i,i-floor(i/3)) = lambda_1+lambda_2;
                phi(i,i-floor(i/3)-1) = lambda_1-lambda_2;
            end    
        end
        phi = phi/2;
    else
        for i = 1:length(x)/2*3
            if mod(i,3) == 1
                phi(i,i-floor(i/3)) = lambda_1;
            elseif mod(i,3) == 2
                phi(i,i-floor(i/3)) = lambda_2;
            end    
        end
    end
    y = phi*x;
elseif test == "shear"
    x1_box = x1;
    x2_box = x1+l1R/2;
    x4_box = x2_box;
    x5_box = x1+l1R;
    x6_box = x5_box+l2R/2;
    x7_box = x5_box+l2R;
    x11_box = x1+l2R;
    x12_box = x1+l2R/2;
    x3_box = x12_box + l1R/2;
    x9_box = x3_box;
    x8_box = x11_box + l1R/2;
    x10_box = x8_box;
    x_box = [x1_box; x2_box; x3_box; x4_box; x5_box; x6_box; x7_box; x8_box; x9_box; x10_box; x11_box; x12_box];
    
    gamma = 0.20;
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

% x and y rigidity constraint matrices (5x12) of (2x2) and (3x3)
L_0 = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
L = [-eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];
if trapezoid % structure with trapezoid-shaped slits
L_0 = [-eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), lambda*eye(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -lambda*eye(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
L = [-eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3);
     eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];
end

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
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-5);
figure
ax = gca;
hold(ax,'on')
set(ax,'YScale','log')
ax.FontSize = 16;
[yOpt, xOpt, Ropt, ax, iterOffset] = minimizationAlgorithm(x, y, Pj, J, R, L, L_0, tol, ax);
titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = false;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");

E = 1;
% check strains
if check_strain
    combos = [1, 2;
              2, 3;
              3, 4;
              4, 5;
              5, 6;
              6, 7;
              7, 8;
              8, 9;
              9, 10;
              10, 11;
              11, 12;
              12, 1;
              3, 12;
              3, 6;
              6, 9;
              9, 12]; % unique node pairings of all edges
    l0 = zeros(length(combos), 1);
    l = zeros(length(combos), 1);
    strains = zeros(length(combos), 1);
    % compute strains for each edge
    for i = 1:length(combos)
        x_a = [xOpt(2*combos(i,1)-1), xOpt(2*combos(i,1))]; % node 1 of the i-th edge
        x_b = [xOpt(2*combos(i,2)-1), xOpt(2*combos(i,2))]; % node 2 of the i-th edge
        l0(i) = norm(x_a - x_b); % reference configuration lengths
        y_a = [yOpt(3*combos(i,1)-2), yOpt(3*combos(i,1)-1), yOpt(3*combos(i,1))];
        y_b = [yOpt(3*combos(i,2)-2), yOpt(3*combos(i,2)-1), yOpt(3*combos(i,2))];
        l(i) = norm(y_a - y_b); % deformed configuration lengths
        strains(i) = (l(i) - l0(i))/l0(i); % resultant strains from the deformation
    end
    [max_val, max_idx] = max(abs(strains));
    disp("Max strain: " + strains(max_idx) + " @ edge (" + combos(max_idx, 1) + ", " + combos(max_idx, 2) + ")")
    disp("Average absolute strain: " + mean(abs(strains)))
    disp("Total absolute strain: " + sum(abs(strains)))
end

% compute geometric stiffness
K = geometricStiffness(xOpt, yOpt, Ropt, Pj, J, L, 2, 2, 2, "translation");
