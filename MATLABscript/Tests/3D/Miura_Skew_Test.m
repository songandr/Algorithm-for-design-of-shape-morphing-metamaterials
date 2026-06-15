% This script is designed as a test case for the MATLAB function 
% minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 06/15/26.
%
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.
close all;
clear all; 
clc;

check_strain = 0;
check_energy = 1;

% Initial x-values
theta_0 = pi/2*0.10;
mu_0 = 1;

l1R = [1; 0];
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

% Define representative parameters @ the centroid of domain in the lambda,
% mu plane, given mu_0 = 1
lambda = (2+cos(theta_0))/3;
mu = (2+cos(theta_0))/3;
theta = acos(cos(theta_0)/(lambda*mu)); % enforce u1 dot v1 = u2 dot v2

figure
hold on
lam = linspace(0,1,500);
mu_curve = cos(theta_0) ./ lam;
idx = mu_curve <= 1;
lam_valid = lam(idx);
mu_valid = mu_curve(idx);
text(0.02,0.98, "$\theta_0 = 0.75(\frac{\pi}{2})$", 'Interpreter', 'latex', "HorizontalAlignment", 'left', 'VerticalAlignment', 'top', 'FontSize', 20)
% Fill region above the curve up to mu=1
fill([lam_valid fliplr(lam_valid)], [mu_valid ones(size(mu_valid))], ...
     [0.4 0.7 0.4], 'EdgeColor','none');
plot(lambda, mu, '.', "MarkerFaceColor", 'k', "MarkerEdgeColor", 'k', "MarkerSize", 6) % centroid point
%text(lambda-0.05, mu, '$(\frac{2+\cos(\theta_0)}{3}, \frac{2+\cos(\theta_0)}{3})$', ...
%    'Interpreter','latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 20)
text(lambda, mu, '$(\lambda_1^*, \lambda_2^*)$', ...
    'Interpreter','latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 20)
plot(lam_valid, mu_valid, 'Color', [0.25 0.25 0.25], 'LineWidth', 2); % boundary curve
text(0.75, cos(theta_0)./0.75, '$\cos(\theta_0)=\lambda_1\lambda_2\cos(\theta)$', 'Interpreter', 'latex','VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 20)
xlim([0 1]); ylim([0 1]);
xticks([0 cos(theta_0) 1]); yticks([cos(theta_0) 1]);
xticklabels({'0', '$\cos(\theta_0)$', '1'}); yticklabels({'$\cos(\theta_0)$', '1'});
xlabel('$\lambda_1$','Interpreter','latex'); ylabel('$\lambda_2$','Interpreter','latex');
ax = gca;
set(ax, 'Layer', 'top');
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';
box on;
axis square;
hold off

l1D = lambda*[1; 0; 0]; % l1D = lambda e_1
l2D = mu*[cos(theta)  -sin(theta)   0; sin(theta) cos(theta)    0;  0   0   1]*[1; 0; 0]; % l2D = mu R(theta) e1

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
L_0 = [-eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2);
     zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), eye(2);
     -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
     zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2);
     zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2);
     eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2)];
    
% y rigidity constraint matrix (7x9) of (3x3)
L = [-eye(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), eye(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), eye(3);
     -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3);
     zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3);
     eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

P1 = [1, 2, 5, 6];
P2 = [2, 3, 4, 5];
P3 = [5, 6, 7, 8];
P4 = [4, 5, 8, 9];

% 3D array containing index set of x coordinates for panel j
Pj = cell(length(J), 1);
for j=1:length(J)
    Pj{j} = eval(sprintf('P%d', j));
end

% Initial R 
for j = 1:length(J)
    R{j} = eye(3); % identity
end

% Initial R 
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
visualizeLatticeVec = true;
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");

E = 1;
count = 0;
max_count = 50;
while E > tol && count < max_count
    % random perturbation post-processing step
    disp("Post-processing...")
    rng(2, "twister");
    std = 10^(-5);
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
    [yOpt_perturbed, xOpt_perturbed, Ropt, ax, iterOffset] = minimizationAlgorithm(x_perturbed, y_perturbed, Pj, J, R, L, L_0, tol, ax, iterOffset);
    vectors = {x_perturbed, y_perturbed, xOpt_perturbed, yOpt_perturbed};
    if E < tol || count == max_count-1
        plot4vectors3D(vectors, titles, visualizeLatticeVec, "Miura");
    end
    xOpt = xOpt_perturbed;
    yOpt = yOpt_perturbed;
    
    % Compute energy associated with solution
    E = 0;

    % Construct necessary cj and rij vectors
    cj = cell(length(J));
    rij = cell(length(J));
    for j = 1:length(J)
        % center of the panel calculation based on new y vector
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
    disp("Energy: " + E)

    count = count + 1;
end

disp("Total post-processing steps: "+count)

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
