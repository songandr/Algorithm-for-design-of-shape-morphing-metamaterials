% This script is designed for running all desired permutations of axial Miura
% tests for the MATLAB function minimizationAlgorithm which is based on the paper:
% 
% "Elastic Energy Approximation and Minimization Algorithm for Foldable
% Meshes"
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 03/19/26.
%
% The initial configuration consists of four panels in the Miura-Ori
% configuration. This test is interested in the folding of these four 
% panels and the resultant energy calculation.

% Initial x-values with lattice vectors e1, e2
s = 0.500025001; % associated with x3-x1 = e1 for gamma = 0.01
gamma = 0.01;
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

lambda = [0.01, 0.25, 0.5, 0.75, 0.9];

y_initial = cell(length(lambda));
xOpt = cell(length(lambda));
yOpt = cell(length(lambda));
x_mathematica = cell(length(lambda));
y_mathematica = cell(length(lambda));
yi_mathematica = cell(length(lambda));
E = cell(length(lambda));

for i=1:length(lambda)
    lambda_1 = lambda(i);
    for j=1:length(lambda)
        disp("i = " + i)
        disp("j = " + j)
        lambda_2 = lambda(j);

        % Construct y from lambda_1, lambda_2
        for k = 1:length(x)/2*3
            if mod(k,3) == 1
                phi(k,k-floor(k/3)) = lambda_1;
            elseif mod(k,3) == 2
                phi(k,k-floor(k/3)) = lambda_2;
            end    
        end
        y = phi*x;
        % Bias y's (y4 y5 y6) in the z-axis
        y(12) = y(12)+0.1;
        y(15) = y(15)+0.1;
        y(18) = y(18)+0.1;

        y_initial{i,j} = y;
        % run Miura axial test for lambda_1, lambda_2
        [xOpt{i,j}, yOpt{i,j}, E{i,j}] = Miura_Axial_Test(x, y);
        x_mathematica{i,j} = makeMathematicaCoords(xOpt{i,j}, "x");
        yi_mathematica{i,j} = makeMathematicaCoords(y_initial{i,j}, "yi");
        y_mathematica{i,j} = makeMathematicaCoords(yOpt{i,j}, "y");
    end
end

xi_mathematica = makeMathematicaCoords(x, "xi");

save('Miura_xi.mat', 'xi_mathematica');
save('Miura_yi.mat', 'yi_mathematica');
save('Miura_x.mat', 'x_mathematica');
save('Miura_y.mat', 'y_mathematica');