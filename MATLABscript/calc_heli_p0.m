function [T1, T2, R1, R2] = calc_heli_p0(q, tau1, z, e)
% 
% returns the operators parameterizing the generators for the helical group
%
% inputs:
% q: an integer satisfying the discreteness condition p*tau1 +
% q*tau2 = 0, p*theta1 + q*theta2 = 2*pi when p = 0.
% tau1: the translation in R3 used to generate the helical group 
% z: a vector in R3 orthogonal to z parameterizing the origin of the isometry
% e: a unit vector in R3 parameterizing the rotation axis
%
% outputs:
% T1: the translation operator in R3 used in g1
% T2: the translation operator in R3 used in g2
% R1: the rotation operator in SO(3) used in g1
% R2: the rotation operator in SO(3) used in g2

I = eye(3,3); % working in R3
ecross = [0,     -e(3),  e(2);
          e(3),  0,      -e(1);
          -e(2), e(1),   0];

% enforcing the discreteness conditions
theta2 = 2*pi/q;
tau2 = 0;

theta1 = 0;

% Rotation matrices about axis e
R1 = I*cos(theta1) + (1-cos(theta1))*(e*e') + sin(theta1)*ecross;
R2 = I*cos(theta2) + (1-cos(theta2))*(e*e') + sin(theta2)*ecross;

% Translation operators
T1 = tau1 * e + (I - R1) * z;
T2 = tau2 * e + (I - R2) * z;

end
