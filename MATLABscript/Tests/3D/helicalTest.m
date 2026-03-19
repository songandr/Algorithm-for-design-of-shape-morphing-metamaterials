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
% Updated Date: 01/20/26.
%
% This test is interested in the folding of panels as a helical origami
% structure and the resultant energy calculation.

% Test Parameters
crease_test = "helical_waterbomb";

e = [0; 1; 0]; % rotation axis: must be unit vector
z = [1; 0; 0]; % origin of isometry: must be orthogonal to e
p = 4;
q = 1;
theta1 = pi/4;

% Parameters for "four_triangles" test
l = 1; % side length of l2R
eta = pi/3; % angle of parallelogram

if crease_test == "helical_four_quads"
    
    % Initial x-values
    l1R = 2*[s; 0];
    l2R = 2*[0; s];
    
    x1 = [2; 0];
    x3 = x1 + l1R;
    x7 = x1 + l2R;
    x9 = x1 + l1R + l2R;
    x2 = (x1+x3)/2;
    x6 = (x1+x7)/2;
    x4 = x6 + l1R;
    x8 = x2 + l2R;
    x5 = (x4+x6)/2;
    % Randomized ICs
    %{
    rng(1, "twister");
    r = normrnd(0, 1/16, [6, 1]); % only perturb eligible DoFs
    x2 = x2 + [r(1); r(2)];
    x8 = x8 + [r(1); r(2)];
    x5 = x5 + [r(3); r(4)];
    x6 = x6 + [r(5); r(6)];
    x4 = x4 + [r(5); r(6)];
    %}
    x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];
    
    % Initial y-values
    
    [T_1, T_2, R_1, R_2] = calc_heli(p, q, theta1, tau1, z, e);
    
    y1 = [r; 0; 0];
    y3 = R_1*y1 + T_1;
    y7 = R_2*y1 + T_2;
    y9 = R_1*y7 + T_1;
    y2 = mean([y1, y3], 2);
    y8 = R_2*y2 + T_2;
    y6 = mean([y1, y7], 2);
    y4 = R_1*y6 + T_1;
    y5 = mean([y1, y9], 2);

    y = [y1; y2; y3; y4; y5; y6; y7; y8; y9];

    % populate the vector numbering all of the panels
    J = 1:4;
    
    % index set for each panel 
    F1 = [1, 2, 5, 6];
    F2 = [2, 3, 4, 5];
    F3 = [5, 6, 7, 8];
    F4 = [4, 5, 8, 9];

    T1 = F1;
    T2 = F2;
    T3 = F3;
    T4 = F4;

    % x rigidity constraint matrix (6x9) of (2x2)
    U = [eye(2),   zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
         eye(2),   zeros(2), -eye(2),  zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
         zeros(2), zeros(2), zeros(2), -eye(2),  zeros(2), eye(2),   zeros(2), zeros(2), zeros(2);
         zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2),   zeros(2), -eye(2);
         eye(2),   zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2),  zeros(2), zeros(2);
         zeros(2), eye(2),   zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), -eye(2),  zeros(2);];
    
    % y rigidity constraint matrix (6x9) of (3x3)
    A = [eye(3),   zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
         R_1,       zeros(3), -eye(3),  zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
         zeros(3), zeros(3), zeros(3), -eye(3),  zeros(3), R_1,       zeros(3), zeros(3), zeros(3);
         zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), R_1,       zeros(3), -eye(3);
         R_2,       zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3),  zeros(3), zeros(3);
         zeros(3), R_2,       zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), -eye(3),  zeros(3);];

elseif crease_test == "helical_four_triangles"
    
    x1 = [2; 0];
    x2 = x1 + [1; 0];
    x4 = x1 + l*[cos(eta); sin(eta)];
    x3 = x4 + [1; 0];
    x0 = mean([x1, x2, x3, x4], 2);
    
    x = [x0; x1; x2; x3; x4];

    tau1 = l;
    [T_1, T_2, R_1, R_2] = calc_heli(p, q, theta1, tau1, z, e);

    y1 = [0; 0; 0];
    y2 = R_1*y1 + T_1;
    y4 = R_2*y1 + T_2;
    y3 = R_2*y2 + T_2;
    y0 = mean([y1, y2, y3, y4], 2);

    y = [y0; y1; y2; y3; y4];

    % populate the vector numbering all of the panels
    J = 1:4;

    % index set for each panel where 1 = x0, 2 = x1, ..., 5 = x4
    F1 = [1, 2, 3];
    F2 = [1, 3, 4];
    F3 = [1, 4, 5];
    F4 = [1, 2, 5];

    T1 = F1;
    T2 = F2;
    T3 = F3;
    T4 = F4;

    % x rigidity constraint matrix (3x5) of (2x2)
    U = [zeros(2), eye(2), zeros(2), zeros(2), zeros(2);
         zeros(2), eye(2),  -eye(2), zeros(2), zeros(2);
         zeros(2), eye(2), zeros(2), zeros(2), -eye(2)];
    
    % y rigidity constraint matrix (3x5) of (3x3)
    A = [zeros(3), eye(3), zeros(3), zeros(3), zeros(3);
         zeros(3), R_1,    -eye(3),   zeros(3), zeros(3);
         zeros(3), R_2,     zeros(3), zeros(3), -eye(3)];

elseif crease_test == "helical_waterbomb"
    
    % relevant parameters
    p = 12; % horizontal chirality index
    q = 3; % vertical chirality index
    %r = 1; % helical radius
    %rng(1, "twister");

    % Initial x-values
    s = 1; % side length
    s2 = s;
    %s = r*sin(pi/p);
    l1R = [s; 0];
    l2R = [0; s2];
    x1 = [0; 0];
    x3 = x1 + l1R;
    x5 = x1 + l2R;
    x2 = (x1+x3)/2;
    %bias1 = normrnd(0, s/64, [2, 1]);
    %x2 = x2 + [0; 1.1*s]; % bias away from midpoint
    x4 = x2 + l2R;
    x6 = x4 - l1R;

    %bias2 = normrnd(0, s/64, [2, 1]);
    %x7 = (x2+x5)/2 + bias2; % bias away from center
    x7 = (x2+x5)/2;
    x = [x1; x2; x3; x4; x5; x6; x7];

    % Initial y-values
    %r = s/(2*sin(pi/p)); % exactly fleshed out helical radius
    r = 0.9*s/(2*sin(pi/p)); % partially folded helical radius

    if q == 0 % degenerate cylinder around y-axis
        e = [0; 1; 0]; % rotation axis: must be unit vector
        z = [0; 0; 0]; % origin of isometry: must be orthogonal to e
        lambda2 = 1;
        l2D = [1, 0;
               0, lambda2;
               0, 0      ] * l2R; % compressed l2R
        [T_1, T_2, R_1, R_2] = calc_heli_q0(p, norm(l2D), z, e);
        theta1 = 2*pi/p;
    elseif p == 0 % degenerate cylinder around x-axis
        e = [1; 0; 0];
        z = [0; 0; 0];
        l1D = [1, 0; 0, 1; 0, 0] * l1R; % compressed l1R
        theta2 = 2*pi/q;
        [T_1, T_2, R_1, R_2] = calc_heli_p0(q, dot(l1D, e), z, e);
    else % helical structure around y-axis
        e = [0; 1; 0];
        z = [0; 0; 0];
        %theta1 = 2*pi/p;
        theta1 = 0.95*norm(x3-x1)/r; % theta1 should be capped and ~l/r
        lambda2 = 0.68; % compression in the e2 direction
        l2D = [1, 0;
               0, lambda2;
               0, 0      ] * l2R; % compressed l2R
        %tau2 = dot(l2D, e);
        tau2 = norm(l2D);
        [T_1, T_2, R_1, R_2] = calc_heli(p, q, theta1, tau2, z, e);
    end

    y1 = r*[1; 0; 0];
    y3 = R_1*y1 + T_1;
    y5 = R_2*y1 + T_2;

    y2 = (y1 + y3)/2;
    y4 = R_2*y2 + T_2;
    y6 = R_1.'*(y4 - T_1);

    y7 = (y2+y5)/2;

    y = [y1; y2; y3; y4; y5; y6; y7];

    % x rigidity constraint matrix (6x7) of (2x2)
    U = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
         -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
         zeros(2), zeros(2), zeros(2), eye(2), zeros(2), -eye(2), zeros(2);
         -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
         zeros(2), -eye(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2)];

    % y rigidity constraint matrix (6x7) of (3x3)
    A = [eye(3), zeros(3), zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3);
         -R_1, zeros(3), eye(3), zeros(3),  zeros(3), zeros(3),  zeros(3);
         zeros(3), zeros(3), zeros(3), eye(3),  zeros(3), -R_1,  zeros(3);
         -R_2, zeros(3), zeros(3), zeros(3),  eye(3), zeros(3),  zeros(3);
         zeros(3), -R_2, zeros(3), eye(3),  zeros(3), zeros(3),  zeros(3)];

    % populate the vector numbering all of the panels
    J = 1:6;

    % index set for each panel where 1 = x0, 2 = x1, ..., 5 = x4
    F1 = [1, 2, 7];
    F2 = [2, 3, 7];
    F3 = [3, 4, 7];
    F4 = [4, 5, 7];
    F5 = [5, 6, 7];
    F6 = [1, 6, 7];

    T1 = F1;
    T2 = F2;
    T3 = F3;
    T4 = F4;
    T5 = F5;
    T6 = F6;

end

% 3D array containing index set of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

Fj = Tj;

% Initial R
R = cell(length(J), 1);
for j = 1:length(J)
    R{j} = eye(3); % identity matrix
end

% initial tolerance for minimization
tol = 10^(-4);

[yOpt, xOpt, Ropt] = minimizationAlgorithmNew(x, y, Fj, Tj, J, R, A, U, tol);

titles = {'Initial X', 'Initial Y', 'Final X', 'Final Y'};
vectors = {x, y, xOpt, yOpt};
visualizeLatticeVec = true;

plot4vectors3D(vectors, titles, visualizeLatticeVec, crease_test);

% check strains
combos = [1, 2;
          2, 3;
          3, 4;
          4, 5;
          5, 6;
          6, 1;
          1, 7;
          2, 7;
          3, 7;
          4, 7;
          5, 7;
          6, 7]; % unique node pairings of all edges
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
