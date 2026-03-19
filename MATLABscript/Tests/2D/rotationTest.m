disp("rotationMin Test")

x1 = [0; 0; 0];
x2 = [1; 0; 0];
x3 = [1; 1; 0];
x4 = [0; 1; 0];
x = [x1; x2; x3; x4];

theta = pi/3;
R = [1  0           0;
     0  cos(theta)  -sin(theta);
     0  sin(theta)  cos(theta)];
disp("theta: ")
disp(theta)
disp("Given rotation matrix: ")
disp(R)

y1 = R*x1;
y2 = R*x2;
y3 = R*x3;
y4 = R*x4;
y = [y1; y2; y3; y4];

c = [0.5; 0.5; 0];
r = [x1-c; x2-c; x3-c; x4-c];

for i = 1:4
       V = r(3*i-2:3*i) + c - y(3*i-2:3*i, 1);
       T = c - r(3*i-2:3*i) - y(3*i-2:3*i);
       Bi{i} =  [0       V(1)    V(2)    V(3);
                T(1)    0       -V(3)   V(2);
                T(2)    V(3)    0       -V(1);
                T(3)    -V(2)   V(1)    0];
end

B = zeros(4);

    for i = 1:4
        B = B + Bi{i}'*Bi{i};
    end

    % rounding the small numeric values to zero
    roundedB = B;
    roundedB(abs(roundedB)<1e-3)=0;

    [V ,D] = eig(roundedB);
    roundedD = D;
    roundedD(abs(roundedD)<1e-3)=0;
    maxEig = max(roundedD, [], "all");
    [~, col] = find(roundedD == maxEig);

    if length(col) > 1
        eigenVal = col(1);
    else
        eigenVal = col;
    end

    pj = V(:, eigenVal);
    v = pj(2:4);
    Rnew = quat2rotm(pj');

Ropt = Rnew;
disp("Optimized rotation matrix: ")
disp(Ropt)