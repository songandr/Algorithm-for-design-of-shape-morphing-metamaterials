% This script is designed for running all desired permutations of axial
% 2D kirigami tests for the MATLAB function minimizationAlgorithm based on 
% the paper:
% 
% "Algorithmic design framework for shape-morphing metamaterials"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 05/03/26.
%
% The initial configuration consists of four panels in the rotating squares
% configuration. This test is interested in the rotation of these four 
% panels with enforced trapezoidal slits and the resultant energy calculation.

% Reminder: set plotCheck = 0 in minimizationAlgorithm

% initialize parameters
b = [-0.25, -0.1, -0.05, 0.05, 0.1, 0.25]; % departure from parallelogram
lambda = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]; % axial shape changes

% Initial x-values
s = 0.5; % side length

l1R = s*[2; 2];
l2R = s*[-2; 2];
l1R = l1R/norm(l1R); % normalize
l2R = l2R/norm(l2R);

x1 = [0; 0]; % P1 square
x2 = [s; 0];
x3 = [s; s];
x12 = [0; s];
x9 = x12 + [0; s];
x4 = [2*s; s]; % P2
x5 = x1 + l1R;
x10 = x2 + l2R; % P4
x11 = [-s; s]; 
x7 = x11 + l1R; % P3
x8 = x4 + l2R;

% y rigidity constraint matrix (5x12) of (3x3)
L = [eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), -eye(3), zeros(3);
     zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3);
     zeros(3), zeros(3), zeros(3), -eye(3), zeros(3), zeros(3), zeros(3), eye(3), zeros(3), zeros(3), zeros(3), zeros(3)];

x_initial = cell(length(b), length(lambda),length(lambda));
y_initial = cell(length(b), length(lambda),length(lambda));
xOpt = cell(length(b), length(lambda),length(lambda));
yOpt = cell(length(b), length(lambda),length(lambda));
E = cell(length(b), length(lambda),length(lambda));

for i=1:length(b)
    b_i = b(i);
    disp("b="+b_i)
    x6 = x3 + [b_i; s]; % construct trapezoid
    x = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];

    % x rigidity constraint matrix (5x12) of (2x2)
    lambda = norm(x6-x9)/norm(x3-x12); % trapezoid parallel side scaling
    L_0 = [eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
         -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2);
         zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), -eye(2), zeros(2);
         zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2);
         zeros(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), zeros(2), zeros(2);
         zeros(2), zeros(2), lambda*eye(2), zeros(2), zeros(2), -eye(2), zeros(2), zeros(2), eye(2), zeros(2), zeros(2), -lambda*eye(2)];

    % perturb initialization to avoid getting stuck
    rng(1, "twister");
    r = normrnd(0, 0.001, [12, 1]); % only 6 DoFs to maintain original lattice vectors post-perturbation
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

    for j=1:length(lambda)
        lambda_1 = lambda(j);
        %disp("lambda1="+lambda_1)
        for k=1:length(lambda)
            lambda_2 = lambda(k);
            %disp("lambda2="+lambda_2)
            % Construct y from lambda_1, lambda_2
            phi = zeros(length(x)/2*3, length(x));
            for l = 1:length(x)/2*3
                if mod(l,3) == 1
                    phi(l,l-floor(l/3)) = lambda_1+lambda_2;
                    phi(l,l-floor(l/3)+1) = lambda_1-lambda_2;
                elseif mod(l,3) == 2
                    phi(l,l-floor(l/3)) = lambda_1+lambda_2;
                    phi(l,l-floor(l/3)-1) = lambda_1-lambda_2;
                end    
            end
            phi = phi/2;
            y = phi*x;
            x_initial{i,j,k} = x;
            y_initial{i,j,k} = y;

            % run Miura axial test for lambda_1, lambda_2
            [xOpt{i,j,k}, yOpt{i,j,k}, E{i,j,k}] = axialKirigamiTest(x,y, L_0, L);

        end
    end
end

% plot results
titles = {'Final X', 'Final Y'};
visualizeLatticeVec = true;
for i=1:8
vectors = {xOpt{1,7,i}, yOpt{1,7,i}};
plot4vectors3D(vectors, titles, visualizeLatticeVec, "Rotating Squares");
end