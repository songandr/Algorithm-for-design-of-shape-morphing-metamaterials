function K = geometricStiffness(xstar, ystar, Rstar, Tj, Fj, J, L, dim, group)
% 
% Based on the paper: "Algorithmic design framework for shape-morphing metamaterials" 
%
% By: Antoine Moats, Yingchao Peng, Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 06/10/26.
%
% Inputs:
% xstar: minimized x coordinate 2-D array (2*n by 1 where n is the number of indices)
% ystar: minimized y coordinate 2-D array (2*n by 1 where n is the number of indices)
% Rstar: minimized rotation matrix for each panel cell array length(J) of 3x3's
% Tj: a column cell array of the set of all x's within each panel 
% Fj: a column cell array of the set of all y's within each panel
% J: the set of all panel labels
% L: matrix of the rigidity constraints (L; chi_1) see Eq. (38)
% dim: 2 for planar rotations or 3 otherwise
% group: "translation" or "helical" for group symmetry category
%
% Outputs:
% K: scalar geometric stiffness where
%   =0: mechanism,
%   >0: bistable

% initialize vectors
n = length(ystar); % number of vertices * 3
lenJ = length(J); % number of panels
xijstar = cell(lenJ, 1); % stores information of x_istar-avg(xstar)_j
yijstar = cell(n/3, lenJ); % stores information of R_jstar*xijstar
% construct chi arrays
[chi, ~] = construct_chi(ystar);
% chi_ij{i,j} stores information of chi_i - avg(chi)_j
chi_ij = cell(n/3, lenJ); % pre-allocate size
% mu_ij{i,j} is an intermediate tensor used for dA computing
mu_ij = cell(3, n);
% dr_jstar is an intermediate tensor used for mu computing
dr_jstar = cell(lenJ, 1);
dA = zeros(n, n);

% construct yijstar and chi_ij
for j = 1:lenJ
    % pos vectors with respect to the center of the panel
    [~, xijstar{j}] = centerOfPanel2D(Tj{j}, xstar);

    num_Pj = length(Fj{j}); % number of vertices in the j-th panel
    chi_avg = (1/num_Pj)*calcMatrixSum_y(chi, Fj{j});
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k

        temp = xijstar{j}(2*i-1:2*i);
        yijstar{k,j} = Rstar{j}*[temp(1); temp(2); 0];

        chi_ij{k,j} = chi{k} - chi_avg; % j dependency encoded in chi_avg
    end
end

% compute mu_ij 2D or 3D
if dim == 2
    e3 = [0; 0; 1]; % out of plane vector

    % compute dummy sums for dr_jstar
    for j = 1:lenJ
        num = zeros(n,1); % numerator for dr_jstar; see Eq. (42, 44)
        den = 0; % denominator for dr_jstar; Eq. (42, 44)
        for i = 1:length(Fj{j}) 
            k = Fj{j}(i); % vertex k
            num = num + chi_ij{k,j}'*cross(e3, yijstar{k,j}); 
            den = den + norm(cross(e3, yijstar{k,j}))^2;
        end
        dr_jstar{j} = num/den;
    end
    
    % compute mu_ij
    for j = 1:lenJ
        for i = 1:length(Fj{j})
            k = Fj{j}(i); % vertex k
            mu_ij{k,j} = chi_ij{k,j} - cross(e3, yijstar{k,j})*dr_jstar{j}';
        end
    end
else
    % compute dummy sums for dr_jstar
    for j = 1:lenJ
        left = zeros(3); % tensor meant to be inverted for dr_jstar; see Eq. (45, 47)
        right = zeros(3, n); % right tensor in dr_jstar; see Eq. (45, 47)
        for i = 1:length(Fj{j})
            k = Fj{j}(i); % vertex k
            yijstarcross = [0                   -yijstar{k,j}(3)    yijstar{k,j}(2);
                            yijstar{k,j}(3)     0                   -yijstar{k,j}(1);
                            -yijstar{k,j}(2)    yijstar{k,j}(1)     0];
            left = left + yijstarcross'*yijstarcross; 
            right = right + yijstarcross'*chi_ij{k,j};
        end
        dr_jstar{j} = left\right;
    end

    % compute mu_ij
    for j = 1:lenJ
        for i = 1:length(Fj{j})
            k = Fj{j}(i); % vertex k
            yijstarcross = [0                   -yijstar{k,j}(3)    yijstar{k,j}(2);
                            yijstar{k,j}(3)     0                   -yijstar{k,j}(1);
                            -yijstar{k,j}(2)    yijstar{k,j}(1)     0];
            mu_ij{k,j} = chi_ij{k,j} - yijstarcross*dr_jstar{j};
        end
    end
end

for j = 1:lenJ
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k
        dA = dA + mu_ij{k,j}'*mu_ij{k,j};
    end
end
tol = 10^(-14);
dA(abs(dA)<tol) = 0; % round off numerics
dA = (dA+dA')/2; % symmetricize

% check for G-H mode special case
N_G = null(L); % null space of (L; chi_1) sized 3I by NG
GH_matrix = N_G'*dA*N_G;
GH_matrix(abs(GH_matrix)<tol) = 0; % round off numerics
GH_matrix = (GH_matrix+GH_matrix')/2; % symmetricize
GH_det = det(GH_matrix);
if GH_det < tol
    K = 0;
    disp("Guest-Hutchinson mode found.")
    return;
end

M_G = (eye(n) - N_G*(GH_matrix\N_G')) * L'; % omit inv(L*L') to use \ operator later
if group == "translation"
    % translation group case
    d1 = L*ystar; d1 = d1(1:3); % first row should be d1
    d2 = L*ystar; d2 = d2(end-5:end-3); % second to last row should be d2
    % last row is always 0 since chi_1 * y = 0 is set
    num_constraints = size(L,1)/3 - 1; % |B_1| + |B_2|; -1 for last row
    
    % minimize group parameters
    T = [d1/norm(d1)   zeros(3,1)  cross(e3,d1)/(sqrt(2)*norm(cross(e3,d1)));
          zeros(3,1)    d2/norm(d2) -cross(e3,d2)/(sqrt(2)*norm(cross(e3,d2)));
          zeros(6,3)]; % see Eq. (55)
    
    G = zeros(3*(num_constraints+1),12); % see Eq. (54); derived from setting g_k = delta d_k
    for k = 1:num_constraints/2 
        G(3*k-2:3*k, 1:3) = eye(3); % assuming 1st half of L holds k=1, 
        G(3*num_constraints/2+3*k-2:3*num_constraints/2+3*k, 4:6) = eye(3); % the 2nd half holds k=2 (before last row)
    end

    B = T'*G'*((L*L')\M_G')*dA*M_G*((L*L')\G)*T; % combining Eq. (53-55) for the final matrix to SVD
    symB = (B+B')/2; % symmetricize
    disp(symB)

    % find the minimum eigenvalue and corresponding eigenvector
    [V, D] = eig(symB); % eigenvectors V and eigenvalues D
    [~,idx] = min(diag(D));
    lambda = V(:,idx);
    lambda = lambda/norm(lambda); % normalize as unit vector
    %disp("(G,0)=")
    %disp(G*T*lambda) % (G,0) = G*T*lambda but also = L*deltay
    % solve for Kstar
    K = dot(lambda,B*lambda);
end
    

if group == "helical"
    % helical group case
end

end