function [yOpt, xOpt, Ropt] = minimizationAlgorithm(x, y, Fj, Tj, J, R, L, L_0, tol)
% 
% Based on the paper: "Elastic Energy Approximation and Minimization
% Algorithm for Foldable Meshes" 
%
% By: Antoine Moats, Niharika Sashidhar, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 03/18/26.
%
% Rigidity Constraints:
% L: matrix of the rigidity constraints that satisfies the equation: L*y = d  
% L_0: matrix of the rigidity constraints that satisfies the equation: L_0*x = c  
%
% Indexing Inputs:
% Tj: a column cell array of the set of all x's within each panel 
% Tj size is 1 by numPanels with cell size vertices in panel by 1;
% x: x coordinate 2-D array (2*n by 1 where n is the number of indices)
% Fj: a column cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth cell)
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
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
%
% Comment out 91-96, 135-141, 224-231 to reduce display clutter
% Comment out 144-158, 235-246 to remove energy plots

% Implement as input eventually
plotCheck = 0;
textOutput = 1;

% initialize cj, rij vectors
n = length(y); % number of vertices * 3
lenJ = length(J);
cj = cell(lenJ);
rij = cell(lenJ);

% construct chi, chihat arrays
[chi, chihat] = construct_chi(y);

% yMap{i,j}, xMap{i,j} stores information of chi_i - avg(chi)_j
yMap = cell(n/3, lenJ); % pre-allocate size
xMap = cell(n/3, lenJ); % pre-allocate size

% initialize A, A_0 arrays for all future vertex minimizations
A = zeros(n, n); % pre-allocating the size of A (3I by 3I)
A_0 = zeros(2*n/3, 2*n/3); % pre-allocating the size of A_0 (2I by 2I)

% construct cj, rij vectors and A, A_0 arrays
for j = 1:lenJ
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel3D(Fj{j}, y);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Tj{j}, x);

    num_Pj = length(Fj{j}); % number of vertices in the j-th panel
    chi_avg = (1/num_Pj)*calcMatrixSum_y(chi, Fj{j});
    chi_0_avg = (1/num_Pj)*calcMatrixSum_x(chihat, Fj{j}); 
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k
        yMap{k, j} = chi{k} - chi_avg;
        xMap{k, j} = chihat{k} - chi_0_avg;
        A = A + yMap{k,j}'*yMap{k,j};
        A_0 = A_0 + xMap{k,j}'*xMap{k,j};
    end
end

% initial invertibility check
%inverseChecks = [];
%inverseChecks(end+1) = inverseCheck(y, Fj, R, A, A_0, L, L_0, yMap, xMap);

% initialize energy E{1}
count = 1;
E{count} = 0;
for j = 1:lenJ
    for i = 1:length(Fj{j})
        k = Fj{j}(i);

        rij_temp = rij{j}(2*i-1:2*i);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);
        
        E{count} = E{count} + norm(y(3*k-2:3*k, 1) - cj{j} - R{j}*[rij1; rij2; 0])^2; 
    end
end 

if textOutput
disp("-------------------------------------")
disp("Iteration number: " + num2str(count));
disp("Initial Energy Value: " + num2str(E{count}));
disp("-------------------------------------")
end

count = count + 1;
E{count} = 0; 

% initial minimizations
RiOpt = minR(x, y, Fj, Tj, J);
yNew = minY(x, y, Tj, J, RiOpt, L, A, yMap);
xNew = minX(x, yNew, Fj, J, RiOpt, L_0, A_0, xMap);

% compute cj, rij
for j = 1:lenJ
    % center of the panel calculation based on new y vector
    [cj{j}, ~] = centerOfPanel3D(Fj{j}, yNew);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Tj{j}, xNew);
end 

% compute the first energy that is to be compared in the while loop (E{2})
for j = 1:lenJ
    for i = 1:length(Fj{j})
        k = Fj{j}(i); 

        rij_temp = rij{j}(2*i-1:2*i);
        rij1 = rij_temp(1);
        rij2 = rij_temp(2);
        
        E{count} = E{count} + norm(yNew(3*k-2:3*k, 1) - cj{j} - RiOpt{j}*[rij1; rij2; 0])^2; 
    end
end 

R = RiOpt;
y = yNew;
x = xNew;

% invertibility check after first minimization
%inverseChecks(end+1) = inverseCheck(y, Fj, R, A, A_0, L, L_0);

err = abs(E{count}-E{count-1});

if textOutput
disp("Iteration number: " + num2str(count));
disp("Previous Energy Value: " + num2str(E{count-1}));
disp("Current Energy Value: " + num2str(E{count}));
disp("Current Error Value: " + num2str(err));
disp("-------------------------------------")
end

% plot energy over iterations
if plotCheck
figure
xlabel("Iteration Number [#]")
ylabel("Energy [L^2]");
%title("Energy Minimization");
set(gca,'YScale','log')
hold on
for i = 1:length(E)
    plot(i-1, E{i}, 'o', 'MarkerFaceColor', [0.10, 0.60, 0.9], 'MarkerEdgeColor', [0.10, 0.60, 0.9]);
    ylim([tol, 10^ceil(log10(max(E{i})))]);
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1])
    drawnow; % ensures that the updated point is plotted
    pause(0.5);
end
end
% initialize while loop parameters
optX = 0; % X needs (not) to be optimized this iteration
optX_count = 1; % counter for # of times x was optimized; see line 66
check = 0; % checker for if x was optimized consecutively
max_x_attempts = 2; % number of times x must be consecutively optimized before breaking loop
optY = 0; % number of times R <-> y has been optimized for the given x
max_y_attempts = 10; % number of times R <-> y can be optimized before forcing an x optimization
force_x = 0; % checker for if x optimization should be forced
finalLoop = 0; % checker for if x is converged and R <-> y must be optimized for final
converged = 0; % checker for if R, y, and x have been optimized (R <-> y is complete once finalLoop = 1)

% optimization loop
while err > tol || converged == 0 % optimize until convergence in R, y, and x
    count = count + 1;
    E{count} = 0;
    
    % only optimize X if R <-> x optimization is converged for the given x
    if optX
        xNew = minX(x, y, Fj, J, R, L_0, A_0, xMap);
        yNew = y;
        optX_count = optX_count+1;
        optY = 0; % reset for a new optimized x
        check = check + 1; % resets to 0 when x doesn't need opt, set to 1 after opt once, set to 2 after consecutive opt
    else % optimize R and y when not converged
        optY = optY + 1;
        R = minR(x, y, Fj, Tj, J);
        yNew = minY(x, y, Tj, J, R, L, A, yMap); 
        xNew = x;
        check = 0;
    end

    if check == max_x_attempts % if X was optimized consecutively max_x_attempts # of times:
        disp("Not converging...")
        break
    end
    
    if optY >= max_y_attempts
        force_x = 1; % force x optimization after max_y_attempts of R <-> y optimization
        if textOutput
        disp("Forcing x optimization...")
        end
    end

    % compute cj and rij vectors
    for j = 1:lenJ
        % center of the panel calculation based on new y vector
        [cj{j}, ~] = centerOfPanel3D(Fj{j}, yNew);

        % pos vectors with respect to the center of the panel
        [~, rij{j}] = centerOfPanel2D(Tj{j}, xNew);
    end

    % compute energy for iteration
    for j = 1:lenJ
        for i = 1:length(Fj{j})
            k = Fj{j}(i); 

            rij_temp = rij{j}(2*i-1:2*i);
            rij1 = rij_temp(1);
            rij2 = rij_temp(2);
        
            E{count} = E{count} + norm(yNew(3*k-2:3*k, 1) - cj{j} - R{j}*[rij1; rij2; 0])^2; 
        end
    end

    err = abs(E{count}-E{count-1});

    % Updating the values to be used at the beginning of the next loop
    y = yNew;
    x = xNew;

    % showing the results of the algorithm in real time
    if textOutput
    disp("Iteration number: " + num2str(count));
    disp("Previous Energy Value: " + num2str(E{count-1}));
    disp("Current Energy Value: " + num2str(E{count}));
    disp("Current Error Value: " + num2str(err));
    disp("-------------------------------------")
    end

    % plotting the results of the algorithm in real time
    if plotCheck
    hold on
    if optX ~= 1
        plot(count-1, E{count}, 'o', 'MarkerFaceColor', [0.10, 0.60, 0.9], 'MarkerEdgeColor', [0.10, 0.60, 0.9]);
    else
        plot(count-1, E{count}, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
    end
    plot(count-1, E{count});
    drawnow; % ensures that the updated point is plotted
    pause(0.1);
    hold off
    end
    % setting boolean parameters for x optimization check
    if (optX == 0 && err < tol) || force_x % x should be optimized next iteration (R <-> y loop converged OR max_y_attempts reached):
        optX = 1;
        force_x = 0; % reset the checker
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
end

Ropt = R;
yOpt = y;
xOpt = x;

disp("x optimization count: " + optX_count);

end
