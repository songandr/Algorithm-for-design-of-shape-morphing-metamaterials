% Planar kirigami designed from recipe; see Fig. 2 from "Mechanism-based
% metamaterials with microstructurally invariant shape-change" Y. Peng et. al

% ICs from Fig. 2a
x1 = [0; 0];
x2 = x1+[1; 0];
x3 = x2+[cos(pi/6); -sin(pi/6)];
x4 = x3+[0;1];
x5 = x2+[0;1];
x6 = x1+[0; 1];
x7 = x6+[cos(pi/3); sin(pi/3)];
x8 = x7+[1;0];
x9 = x4+[cos(pi/3); sin(pi/3)];
x = [x1; x2; x3; x4; x5; x6; x7; x8; x9];
l1R_0 = x3-x1;
l2R_0 = x7-x1;

% helper panels
x_right = x;
x_up = x;
for i=1:9
    x_right(2*i-1:2*i) = x(2*i-1:2*i) + l1R_0;
    x_up(2*i-1:2*i) = x(2*i-1:2*i) + l2R_0;
end

% populate the vector numbering all of the panels
J = [1, 2, 3, 4];

% index set for panels 1~4 
T1 = [1, 2, 5, 6];
T2 = [2, 3, 4, 5];
T3 = [5, 6, 7];
T4 = [4, 5, 7, 8, 9];

% cell array containing labeling of x coordinates for panel j
Tj = cell(length(J), 1);
for j=1:length(J)
    Tj{j} = eval(sprintf('T%d', j));
end

% initialize cj vectors
cj = cell(length(J));
cj_right = cell(length(J));
cj_up = cell(length(J));
for j=1:length(J)
    % center of the panel calculation based on initial x vector
    [cj{j}, ~] = centerOfPanel2D(Tj{j}, x);
    [cj_right{j}, ~] = centerOfPanel2D(Tj{j}, x_right);
    [cj_up{j}, ~] = centerOfPanel2D(Tj{j}, x_up);

end 

% transform centroids by tensor D and necessary panel translations
x_new = zeros(2*sum(cellfun(@numel, Tj)),1); % number of vertices for new structure is number of vertices in Tj; append new panels later
x_new_right = zeros(2*sum(cellfun(@numel, Tj)),1);
x_new_up = zeros(2*sum(cellfun(@numel, Tj)),1);
cj_new = cell(length(J));
rj_new = cell(length(J));
cj_new_right = cell(length(J));
rj_new_right = cell(length(J));
cj_new_up = cell(length(J));
rj_new_up = cell(length(J));

d1 = 0; % d1 = d2 = 0 -> conformal case
d2 = 0; % must be in the region of validity (Fig. S4)
D = 1/(1+4*d1*d2) * [2       4*d2;
                     -4*d1   2   ]; % diagonal shape tensor from Eq. S9
k = 1;
for j=1:length(J)
    cj_new{j} = D*cj{j};
    rj_new{j} = cj_new{j} - cj{j};
    cj_new_right{j} = D*cj_right{j};
    rj_new_right{j} = cj_new_right{j} - cj_right{j};
    cj_new_up{j} = D*cj_up{j};
    rj_new_up{j} = cj_new_up{j} - cj_up{j};
    for i=1:length(Tj{j})
        x_i = [x(2*i-1); x(2*i)];
        x_new(k:k+1) = x_i + rj_new{j};

        x_i_right = [x_right(2*i-1); x_right(2*i)];
        x_new_right(k:k+1) = x_i_right + rj_new_right{j};
        x_i_up = [x_up(2*i-1); x_up(2*i)];
        x_new_up(k:k+1) = x_i_up + rj_new_up{j};

        k = k+2;
    end
end

% add connecting panels between spaced out panels
J_added = 0; % number of new panels = number of shared edges
for j=1:length(J)-1 % 1st panel for comparison
    for i=1:length(Tj{j}) % iterate over vertices
        vertex = Tj{j}(i); % vertex of interest
        for k=j+1:length(J) % other panels to compare with
            if ismember(vertex, Tj{k}) % one shared vertex
                other_vertices = Tj{j}(Tj{j} ~= vertex);
                for l=1:length(other_vertices)
                    vertex2 = other_vertices(l);
                    if ismember(vertex2, Tj{k}) % two shared vertices
                        J_added = J_added + 1;
                    end
                end
            end
        end
    end
end
J_added = J_added/2; % adjust for double count

% relabel transformed panels and vertices
Tj_new = [Tj; cell(J_added, 1)];
count = 1;
for j=1:length(Tj)
    for i=1:length(Tj{j})
        if count > sum(cellfun(@numel, Tj))
            break
        end
        Tj_new{j}(i) = count;
        count = count+1;
    end
end

% add remaining connector vertices
l1R = x_new_right(1:2) - x_new(1:2);
l2R = x_new_up(1:2) - x_new(1:2);
x17 = x_new(15*2-1:15*2) - l2R;
x18 = x_new(14*2-1:14*2) - l2R;
x19 = x_new(11*2-1:11*2) + l1R - l2R;
x20 = x17 + l1R - (x_new(3:4) - x_new(1:2)); % parallelogram
x21 = x_new(1:2) + l1R;
x22 = x_new(4*2-1:4*2) + l1R;
x23 = x_new(9*2-1:9*2) + l1R;

x_new = [x_new; x17; x18; x19; x20; x21; x22; x23];

% add connector panels
Tj_new{5} = [2, 5, 17];
Tj_new{6} = [6, 18, 19, 20, 21];
Tj_new{7} = [3, 8, 10, 17];
Tj_new{8} = [7, 13, 22, 23];

vector = {x_new};
visualizeLatticeVec = 0;
plot1vector(vector,'', 'Initial X', visualizeLatticeVec, "Planar Kirigami");
