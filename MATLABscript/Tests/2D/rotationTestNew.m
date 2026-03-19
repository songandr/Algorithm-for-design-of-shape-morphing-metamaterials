disp("rotationMin_new Test")

x1 = [0; 0; 0];
x2 = [1; 0; 0];
x3 = [1; 1; 0];
x4 = [0; 1; 0];
x = [x1; x2; x3; x4];

theta = pi/4;
R = [cos(theta)  0  sin(theta);
     0  1  0;
     -sin(theta)    0  cos(theta)];
disp("theta: ")
disp(theta)
disp("Given rotation matrix: ")
disp(R)

lambda = [1 0 0;
          0 1 0;
          0 0 2];

y1 = R*lambda*x1;
y2 = R*lambda*x2;
y3 = R*lambda*x3;
y4 = R*lambda*x4;
y = [y1; y2; y3; y4];

c = [0.5; 0.5; 0];
rkj = [x1-c; x2-c; x3-c; x4-c];

Ml = zeros(4, 4); % Ml = sum of sym(M_kl) over k (1 panel)

for i = 1:4
    ckl = y(3*i-2:3*i) - c; 
    ckl_cross = [0           -ckl(3)    ckl(2);
                 ckl(3)      0          -ckl(1);
                 -ckl(2)     ckl(1)     0];

    r = rkj(3*i-2:3*i);

    rkl_cross = [0       -r(3) r(2);
                 r(3)    0     -r(1);
                 -r(2)   r(1)  0];

    Mkl11 = dot(r, ckl);
    Mkl12 = transpose(cross(r, ckl));
    Mkl21 = cross(r, ckl);
    Mkl22 = ckl*r.' - transpose(rkl_cross)*(ckl_cross);
    Mkl22 = Mkl22;

    Mkl{i} = [Mkl11    Mkl12(1)   Mkl12(2)   Mkl12(3);
              Mkl21(1) Mkl22(1,1) Mkl22(1,2) Mkl22(1,3);
              Mkl21(2) Mkl22(2,1) Mkl22(2,2) Mkl22(2,3);
              Mkl21(3) Mkl22(3,1) Mkl22(3,2) Mkl22(3,3)];

    Mkl_sym{i} = (Mkl{i}+Mkl{i}.')/2;
    Ml = Ml + Mkl_sym{i};
end

roundedMl = Ml;
roundedMl(abs(roundedMl)<1e-4) = 0;

[V, D] = eig(roundedMl); % eigenvectors V and eigenvalues D
    
% rounding the small numeric values to zero
roundedD = D;
roundedD(abs(roundedD)<1e-3)=0;

% find the maximum eigenvalue and corresponding eigenvector @ col
maxEig = max(roundedD, [], "all");
[~, col] = find(roundedD == maxEig);
    
if length(col) > 1
    eigenVal = col(1);
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
R = v*v.' + qr^2*eye(3) + 2*qr*v_cross + v_cross*v_cross;

disp("Optimized rotation matrix: ")
disp(R)