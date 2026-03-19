function [xOpt, yOpt, E] = Miura_Axial_Test(x, y)
% 
% runs 1 axial shape change test on a Miura topology with N vertices
%
% inputs:
% x: a 2N sized column array holding reference configuration vertices
% y: a 3N sized column array holding deformed configuration vertices
%
% outputs:
% xOpt: energy minimizing reference configuration vertices (length 2N)
% yOpt: energy minimizing deformed configuration vertices (length 3N)

% x rigidity constraint matrix (7x9) of (2x2)
L_0 = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
    -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);];
    
% y rigidity constraint matrix (7x9) of (3x3)
L = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
    -eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);];

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
for j = 1:length(J)
    R{j} = eye(3); % identity
end

% tolerance for minimization
tol = 10^(-5);

[yOpt, xOpt, ~] = minimizationAlgorithm(x, y, Fj, Tj, J, R, L, L_0, tol);

% random perturbation post-processing step
E = 1; % initialize energy
while E > 10^(-4)

r = normrnd(0, 0.01, [length(x), 1]);
r_y = normrnd(0, 0.01, [length(y), 1]);
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
[yOpt, xOpt, Ropt] = minimizationAlgorithm(x_perturbed, y_perturbed, Fj, Tj, J, R, L, L_0, tol);

% Compute energy associated with solution
E = 0;

% Construct necessary cj and rij vectors
cj = cell(length(J));
rij = cell(length(J));
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel3D(Fj{j}, yOpt);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Tj{j}, xOpt);
end

for j = 1:length(J)
    for i = 1:length(Fj{j})
        k = Fj{j}(i);

        rij_temp = rij{j}(2*i-1:2*i);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);
        
        E = E + norm(yOpt(3*k-2:3*k, 1) - cj{j} - Ropt{j}*[rij1; rij2; 0])^2;
    end
end 

end
disp("Energy: " + E)
end