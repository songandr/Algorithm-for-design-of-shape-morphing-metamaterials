function Ropt = minR(x, y, Pj, J)
% 
% uses quaternions to convert the minimization into an eigenvalue problem, 
% the iterative rotation minimization determines the rotation matrix for each panel
% that yields a convergent, minimized energy based on a configuration of
% given X and Y values
%
% inputs:
% x: x coordinate 2-D array (2*n by 1 where n is the number of indices)
% y: y coordinate 2-D array (3*n by 1 where n is the number of indices)
% Pj: a column cell array of the set of all y's within each panel (the jth panel
% corresponds to the jth cell)
% J: the set of all panels
%
% outputs:
% Ropt: a cell array of the rotation matrices that minimize the elastic energy 

% setting the values of the initial cj and rij vectors

cj = cell(length(J));
rij = cell(length(J));
for j = 1:length(J)
    % center of the panel calculation based on initial y vector
    [cj{j}, ~] = centerOfPanel3D(Pj{j}, y);

    % pos vectors with respect to the center of the panel
    [~, rij{j}] = centerOfPanel2D(Pj{j}, x);
end 

% initialize rotation matrix
Ml = zeros(4, 4, length(J)); % Ml = sum of sym(M_kl) over k in panel j

% construct the matrix to optimize over
for j = 1:length(J)

    for i = 1:length(Pj{j})

       k = Pj{j}(i);
       ckl = y(3*k-2:3*k, 1) - cj{j}; 
       ckl_cross = [0           -ckl(3)    ckl(2);
                    ckl(3)      0          -ckl(1);
                    -ckl(2)     ckl(1)     0];
       
       r_temp = rij{j}(2*i-1:2*i);
       r = [r_temp(1); r_temp(2); 0]; % convert to R3 for SO3 rotation compatibility

       rkl_cross = [0       -r(3) r(2);
                    r(3)    0     -r(1);
                    -r(2)   r(1)  0];

       Mkl11 = dot(r, ckl);
       Mkl12 = transpose(cross(r, ckl));
       Mkl21 = cross(r, ckl);
       Mkl22 = ckl*r.' - transpose(rkl_cross)*(ckl_cross);

       Mkl = [Mkl11    Mkl12(1)   Mkl12(2)   Mkl12(3);
              Mkl21(1) Mkl22(1,1) Mkl22(1,2) Mkl22(1,3);
              Mkl21(2) Mkl22(2,1) Mkl22(2,2) Mkl22(2,3);
              Mkl21(3) Mkl22(3,1) Mkl22(3,2) Mkl22(3,3)];

       Mkl_sym = (Mkl+Mkl.')/2;
       Ml(:,:,j) = Ml(:,:,j) + Mkl_sym;
    end
end

% optimize the quaternion and construct the rotation matrix

for j = 1:length(J)

    % rounding small numeric values to zero
    roundedMl = Ml(:,:,j);
    roundedMl(abs(roundedMl)<1e-10) = 0;

    [V, D] = eig(roundedMl); % eigenvectors V and eigenvalues D
    
    % rounding small numeric values to zero
    roundedD = D;
    roundedD(abs(roundedD)<1e-10)=0;

    % find the maximum eigenvalue and corresponding eigenvector @ col
    maxEig = max(roundedD, [], "all");
    [~, col] = find(roundedD == maxEig);
    
    if length(col) > 1 % eigenvalue multiplicity check
        error('Error: multiple rotation eigenvalues found...');
    else
        eigenVal = col;
    end

    q = V(:, eigenVal); % eigenvector associated with greatest eigenvalue
    q = q/norm(q); % normalized

    % unpack quaternion form into a rotation tensor
    qr = q(1);
    v = q(2:4);
    v_cross = [0       -v(3) v(2);
               v(3)    0     -v(1);
               -v(2)   v(1)  0];
    R{j} = v*v.' + qr^2*eye(3) + 2*qr*v_cross + v_cross*v_cross;
end

Ropt = R;
end 