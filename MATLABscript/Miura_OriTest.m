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
% Updated Date: 06/12/26.
%
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.

test = "axial"; % [reference, axial, shear]. Shear tests must have non-orthogonal reference lattice vectors.
check_strain = 0;
check_energy = 1;

% Initial x-values
%s = 0.500025001;  % associated with x3-x1 = e1 for gamma = 0.01
%s = 0.5025104592;  % associated with x3-x1 = e1 for gamma = 0.1
%gamma = 0.1;
s = 0.5025104592;
gamma = 0.1;
x1 = [0; 0];
x2 = s*[cos(gamma); sin(gamma)];
x3 = s*[2*cos(gamma); 0];
x4 = s*[2*cos(gamma); 0] + [0; 0.5];
x5 = s*[cos(gamma); sin(gamma)] + [0; 0.5];
x6 = [0; 0.5];
x7 = [0; 1];
x8 = x2 + [0; 1];
x9 = x3 + [0; 1];
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];

phi = zeros(length(x)/2*3, length(x));
if test == "reference"
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = 1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = 1;
        end
    end
    
elseif test == "axial"
    lambda_1 = 0.8; % axial deformation
    lambda_2 = 0.8;
    for i = 1:length(x)/2*3
        if mod(i,3) == 1
            phi(i,i-floor(i/3)) = lambda_1;
        elseif mod(i,3) == 2
            phi(i,i-floor(i/3)) = lambda_2;
        end    
    end

elseif test == "shear"
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
end

y = phi*x;

% Bias y's (y4 y5 y6) in the z-axis
y(12) = y(12)+0.1;
y(15) = y(15)+0.1;
y(18) = y(18)+0.1;

% x rigidity constraint matrix (7x9) of (2x2)
L_0 = [-eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
    
% y rigidity constraint matrix (7x9) of (3x3)
L = [-eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

T1 = [1, 2, 5, 6];
T2 = [2, 3, 4, 5];
T3 = [5, 6, 7, 8];
T4 = [4, 5, 8, 9];

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R 
R = cell(1, length(J));
for j = 1:length(J)
    R{j} = eye(3); % identity
end

% initial tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, ~] = minimizationAlgorithm(x, y, Fj, Tj, J, R, L, L_0, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;
%{
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");
%}
E = 1;
count = 0;
%while E > 2*10^(-5) && count < 2
while count < 2
% random perturbation post-processing step
disp("Post-processing...")
rng(2, "twister");
std = 1/32;
if count == 1
    std = 1/128;
end
r = normrnd(0, std, [length(x), 1]);
r_y = normrnd(0, std, [length(y), 1]);
% only perturb eligible DoFs
r(1:2) = 0; % fixed first node
r(17:18) = 0; % top right corner
r(5:6) = 0; % bottom right corner
r(13:14) = 0; % top left corner
r_y(1:3) = 0; % fixed first node
r_y(end-2:end) = 0; % top right corner
r_y(7:9) = 0; % bottom right corner
r_y(end-8:end-6) = 0; % top left corner

x_perturbed = xOpt + r;
y_perturbed = yOpt + r_y;

% enforce original constraints post-perturbation
x_perturbed(7) = xOpt(7) + r(11); % shift right side
x_perturbed(8) = xOpt(8) + r(12);
x_perturbed(15) = xOpt(15) + r(3); % shift top side
x_perturbed(16) = xOpt(16) + r(4);

y_perturbed(10) = yOpt(10) + r_y(16); % shift right side
y_perturbed(11) = yOpt(11) + r_y(17);
y_perturbed(12) = yOpt(12) + r_y(18);
y_perturbed(end-5) = yOpt(end-5) + r_y(4); % shift top side
y_perturbed(end-4) = yOpt(end-4) + r_y(5);
y_perturbed(end-3) = yOpt(end-3) + r_y(6);

disp("Starting perturbation minimization...")
[yOpt_perturbed, xOpt_perturbed, Ropt] = minimizationAlgorithm(x_perturbed, y_perturbed, Fj, Tj, J, R, L, L_0, tol);
vectors = {x_perturbed, y_perturbed, xOpt_perturbed, yOpt_perturbed};
if count == 1
    plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");
end
xOpt = xOpt_perturbed;
yOpt = yOpt_perturbed;
%{
xi_mathematica = makeMathematicaCoords(x, "xi");
yi_mathematica = makeMathematicaCoords(y, "yi");
x_mathematica = makeMathematicaCoords(xOpt_perturbed, "x");
y_mathematica = makeMathematicaCoords(yOpt_perturbed, "y");
%}
if check_energy
    % compute cj, rij
    cj = cell(length(J));
    rij = cell(length(J));
    for j = 1:length(J)
        % center of the panel calculation based on new y vector
        [cj{j}, ~] = centerOfPanel3D(Fj{j}, yOpt_perturbed);
    
        % pos vectors with respect to the center of the panel
        [~, rij{j}] = centerOfPanel2D(Tj{j}, xOpt_perturbed);
    end 
    
    % compute energy
    E = 0;
    for j = 1:length(J)
        for i = 1:length(Fj{j})
            k = Fj{j}(i);
            
            rij_temp = rij{j}(2*i-1:2*i);
            rij1 = rij_temp(1);
            rij2 = rij_temp(2);
            
            E = E + norm(yOpt_perturbed(3*k-2:3*k, 1) - cj{j} - Ropt{j}*[rij1; rij2; 0])^2;
        end
    end 
    disp("Energy: " + E)
end
count = count + 1;
end

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
              1, 6;
              2, 5;
              5, 8;
              4, 9]; % unique node pairings of all edges
    l0 = zeros(length(combos), 1);
    l = zeros(length(combos), 1);
    strains = zeros(length(combos), 1);
    % compute strains for each edge
    for i = 1:length(combos)
        x_a = [xOpt_perturbed(2*combos(i,1)-1), xOpt_perturbed(2*combos(i,1))]; % node 1 of the i-th edge
        x_b = [xOpt_perturbed(2*combos(i,2)-1), xOpt_perturbed(2*combos(i,2))]; % node 2 of the i-th edge
        l0(i) = norm(x_a - x_b); % reference configuration lengths
        y_a = [yOpt_perturbed(3*combos(i,1)-2), yOpt_perturbed(3*combos(i,1)-1), yOpt(3*combos(i,1))];
        y_b = [yOpt_perturbed(3*combos(i,2)-2), yOpt_perturbed(3*combos(i,2)-1), yOpt(3*combos(i,2))];
        l(i) = norm(y_a - y_b); % deformed configuration lengths
        strains(i) = (l(i) - l0(i))/l0(i); % resultant strains from the deformation
    end
    [max_val, max_idx] = max(abs(strains));
    disp("Max strain: " + strains(max_idx) + " @ edge (" + combos(max_idx, 1) + ", " + combos(max_idx, 2) + ")")
    disp("Average absolute strain: " + mean(abs(strains)))
    disp("Total absolute strain: " + sum(abs(strains)))
end

% compute geometric stiffness
K = geometricStiffness(xOpt_perturbed, yOpt_perturbed, Ropt, Tj, Fj, J, L, 3, 2, 3, "translation");

