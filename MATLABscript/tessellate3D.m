function x_tessellated = tessellate3D(x, l1, l2, n)

% Function to tessellate a set of vertices in R3 by Bravais lattice vectors
%
% By: Andrew Song
% Under the Supervision of Dr. Paul Plucinsky
% Viterbi School of Engineering, Unversity of Southern California 
%
% Updated Date: 03/26/25
%
% Inputs:
% x: coordinates of vertices, 3N by 1 vector where N is # of vertices
% l1: horizontal Bravais lattice vector [l1_x l1_y 0]
% l2: vertical Bravais lattice vector [l2_x l2_y 0]
% n: # of tessellations (if n=0 -> tessellate3D is identity)

% convert Bravais lattice vector to 1 by 2 to vectorize operations:
l1 = [l1(1) l1(2) l1(3)];
l2 = [l2(1) l2(2) l1(3)];

% convert x to N by 2 array
N = length(x)/3; % number of vertices
x_new = zeros(N, 3);
for i = 1:N
    x_new(i,1) = x(3*i-2, 1);
    x_new(i,2) = x(3*i-1, 1);
    x_new(i,3) = x(3*i, 1);
end

% initialize output array (to be reorganized to vector)
x_tessellated = zeros((n+1)^2*N,3);
x_tessellated(1:N, :) = x_new;

% tessellate n times
for i = 1:N % vertices
    for j = 1:n % tessellation loops (1st is a 2x2 tessellation)
        for k = 1:j % total of j horizontal, vertical tesselations
            if j ~= k % tessellate from non-corner cells
                % counting index in terms of j and k
                current_index_1 = j^2 + (k-1)*2+1; % applicable for horizontal tessellations
                current_index_2 = j^2 + 2*k; % applicable for vertical tessellations
                prev_index_1 = (j-1)^2 + (k-1)*2+1; % index being tessellated horizontally
                prev_index_2 = (j-1)^2 + 2*k; % index being tessellated horizontally
                
                % tessellate horizontally
                x_tessellated((current_index_1-1)*N+i, :) = x_tessellated((prev_index_1-1)*N+i, :) + l1;

                % tessellate vertically
                x_tessellated((current_index_2-1)*N+i, :) = x_tessellated((prev_index_2-1)*N+i, :) + l2;

            else % tessellate from last corner cell
                % tessellate horizontally
                x_tessellated(N*((j+1)^2 - 3)+i, :) = x_tessellated(N*(j^2-1)+i,:) + l1;
        
                % tessellate vertically
                x_tessellated(N*((j+1)^2 - 2)+i, :) = x_tessellated(N*(j^2-1)+i,:) + l2;

            end
        end

        % tessellate diagonally from previous corner
        x_tessellated(((j+1)^2-1)*N+i, :) = x_tessellated(N*(j^2-1)+i, :) + l1 + l2;

    end
end

% remove redundant vertices
x_tessellated = unique(x_tessellated, 'rows');

%{
% reorganize x_tessellated to original column vector structure
x_final = zeros(length(x_tessellated)*2, 1);
for i = 1:length(x_tessellated)
    x_final(2*i-1, 1) = x_tessellated(i, 1);
    x_final(2*i, 1) = x_tessellated(i, 2);
end

x_tessellated = x_final; % output
%}
end