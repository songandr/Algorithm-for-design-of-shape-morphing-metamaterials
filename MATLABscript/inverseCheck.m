function inverseCheck = inverseCheck(y, Fj, R, A, A_0, L, L_0, yMap, xMap)
% 
% checks the invertibility for a simultaneous optimization of x and y given
% a set of panel rotations
%
% Inputs
% y: y coordinate 2-D array (3*n by 1 where n is the number of vertices)
% Fj: a column cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth row)
% R: a cell array of all of the rotation matrices for each panel
% A, A_0: precomputed
% L, L_0: matrix of the rigidity constraints that satisfies the equation:
%         L*y = d, L_0*x =c
% Outputs
% inverseCheck: boolean on whether the required matrix is invertible

% Initialize relevant parameters
n = length(y); % number of vertices * 3
lenJ = length(Fj); % number of panels

A_cross = zeros(n, 2*n/3); % pre-allocating the size of A_cross (3I by 2I)
for j = 1:lenJ
    for i = 1:length(Fj{j})
        k = Fj{j}(i); % vertex k
        A_cross = A_cross - yMap{k,j}'*R{j}*xMap{k,j};
    end
end

N_G = null(L);
N_T = null(L_0);

sym_1 = 0.5*(N_G'*A*N_G + (N_G'*A*N_G)');
sym_2 = 0.5*(N_T'*A_0*N_T + (N_T'*A_0*N_T)');

D = [sym_1                  N_G'*A_cross*N_T;
     (N_G'*A_cross*N_T)'    sym_2];

%D = eye(size(N_T,2)) - (N_T'*A_0*N_T)\(N_T'*A_cross'*N_G) * ((N_G'*A*N_G)\(N_T'*A_cross'*N_G)');
%D_sym = 0.5*(D+D');
%D_sym = chop(5, D_sym);


if rcond(D) < 1e10 % singular
    inverseCheck = 0;
else % invertible
    inverseCheck = 1;
end



