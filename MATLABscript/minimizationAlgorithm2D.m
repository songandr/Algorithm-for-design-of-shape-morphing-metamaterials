function [yOpt, xOpt, Ropt] = minimizationAlgorithm2D(x, y, Fj, Tj, J, R, A, U, tol)
% 
% Based on the paper: "Elastic Energy Approximation and Minimization
% Algorithm for Foldable Meshes" 
%
% By: Antoine Moats, Niharika Sashidhar, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 03/26/25
%
% Rigidity Constraints:
% A: matrix of the rigidity constraints that satisfies the equation: Ay = e  
% U: U is a matrix that contains the rigidity constraints such that Ux = h
%
% Indexing Inputs:
% Tj: a 3D array of the set of all x's within each panel 
% Tj size is 1 by verticesPerPanel  by numPanels;
% x: x coordinate 2-D array (2*n by 1 where n is the number of indices)
% Fj: a 2D array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% y: y coordinate 2-D array (2*n by 1 where n is the number of indices)
% J: the set of all panels
% R: a cell array of all of the initial rotation matrix for each panel
% tol: tolerance for the algorithm minimization
%
% Outputs:
% yMin: y coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% xMin: x coordinate 2-D array that minimizes the elastic energy based on
% given rigidity constraints  
% Ropt: array of rotation matrices for minimizes the elastic energy  

% initialize cj and rij vectors
rij = zeros(2*length(Tj), 1, length(J));
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel2D(Fj(:, :, j), y);

    % pos vectors with respect to the center of the panel
    [~, rij(:, :, j)] = centerOfPanel2D(Tj(:, :, j), x);
end 

% initialize energy E{1}
count = 1;
E{count} = 0;
for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
        k = Fj(:, i, j);

        rij_temp = rij(2*i-1:2*i, 1, j);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);

        E{count} = E{count} + norm(y(2*k-1:2*k, 1) - cj{j} - R{j}*[rij1; rij2])^2;
    end
end 

count = count + 1;
E{count} = 0; 

% initial minimizations
RiOpt = rotationMin_2D(x, y, Fj, Tj, J);
yNew = minY_V2_2D(x, y, Fj, Tj, J, RiOpt, A);
xNew = minX_V2(x, yNew, Fj, Tj, J, RiOpt, U);

% compute cj, rij
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel2D(Fj(:, :, j), yNew);

    % pos vectors with respect to the center of the panel
    [~, rij(:, :, j)] = centerOfPanel2D(Tj(:, :, j), xNew);
end 

% compute the first energy that is to be compared in the while loop (E{2})
for j = 1:length(J)
    for i = 1:length(Fj(:, :, j))
        k = Fj(:, i, j); 

        rij_temp = rij(2*i-1:2*i, 1, j);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);

        E{count} = E{count} + norm(yNew(2*k-1:2*k, 1) - cj{j} - RiOpt{j}*[rij1; rij2])^2;
    end
end 

R = RiOpt;
y = yNew;
x = xNew;

err = abs(E{count}-E{count-1});

% plot energy over iterations
figure
xlabel("Iteration Number [#]")
ylabel("Energy [L^2]");
title("Energy Minimization");

hold on
for i = 1:length(E)
    scatter(i-1, E{i}, 'filled', 'MarkerFaceColor', [0.10, 0.60, 0.9]);
    drawnow; % ensures that the updated point is plotted
    pause(0.5);
end

% initialize while loop parameters
optX = 0; % X needs (not) to be optimized this iteration
optX_count = 1; % counter for # of times x was optimized
check = 0; % checker for if x was optimized consecutively
max_attempts = 2; % number of times x must be consecutively optimized before breaking loop
finalLoop = 0; % checker for if x is converged and R <-> y must be optimized for final
converged = 0; % checker for if R, y, and x have been optimized (R <-> y is complete once finalLoop = 1)

% optimization loop
while err > tol || converged == 0 % optimize until convergence in R, y, and x
    count = count + 1;
    E{count} = 0;
    
    % only optimize X if R <-> x optimization is converged for the given x
    if optX
        xNew = minX_V2(x, y, Fj, Tj, J, R, U);
        yNew = y;
        optX_count = optX_count+1;
        check = check + 1; % resets to 0 when x doesn't need opt, set to 1 after opt once, set to 2 after consecutive opt
    else % optimize R and y when not converged
        R = rotationMin_2D(x, y, Fj, Tj, J);
        yNew = minY_V2_2D(x, y, Fj, Tj, J, R, A); 
        xNew = x;
        check = 0;
    end

    if check == max_attempts % if X was optimized consecutively max_attempts # of times:
        disp("Not converging...")
        break
    end

    % compute cj and rij vectors
    for j = 1:length(J)
        % center of the panel calculation based on y vector
        [cj{j}, ~] = centerOfPanel2D(Fj(:, :, j), yNew);

        % pos vectors with respect to the center of the panel
        [~, rij(:, :, j)] = centerOfPanel2D(Tj(:, :, j), xNew);
    end

    % compute energy for iteration
    for j = 1:length(J)
        for i = 1:length(Fj(:, :, j))
            k = Fj(:, i, j); 

            rij_temp = rij(2*i-1:2*i, 1, j);
            rij1 = rij_temp(1);
            rij2 = rij_temp(2);
        
            E{count} = E{count} + norm(yNew(2*k-1:2*k, 1) - cj{j} - R{j}*[rij1; rij2])^2; 
        end
    end

    err = abs(E{count}-E{count-1});

    % Updating the values to be used at the beginning of the next loop
    y = yNew;
    x = xNew;

    % showing the results of the algorithm in real time
    disp("-------------------------------------")
    disp("Iteration number: " + num2str(count));
    disp("Total Energy Calculations: " + num2str(count));
    disp("-------------------------------------")
    disp("Previous Energy Value: " + num2str(E{count-1}));
    disp("Current Energy Value: " + num2str(E{count}));
    disp("Current Error Value: " + num2str(err));

    % setting boolean parameters for x optimization check
    if optX == 0 && err < tol % x should be optimized next iteration (R <-> y loop converged):
        optX = 1;
        if finalLoop % energy has converged for R, y, and x minimizations
            disp("Converged!")
            converged = 1;
        end
    elseif optX == 1 && err < tol % R <-> y should be optimized one last time (x converged):
        optX = 0;
        finalLoop = 1;
    else % continue optimizing R <-> y
        optX = 0;
    end
    
    % plotting the results of the algorithm in real time
    hold on
    scatter(count-1, E{count}, 'filled', 'MarkerFaceColor', [0.10, 0.60, 0.9]);
    plot(count-1, E{count});
    drawnow; % ensures that the updated point is plotted
    pause(0.1); %pausing for 1/2 of a second
    hold off
end

Ropt = R;
yOpt = y;
xOpt = x;

disp("x optimization count: " + optX_count);

end
