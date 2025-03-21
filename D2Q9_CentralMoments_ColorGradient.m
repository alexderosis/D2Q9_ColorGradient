%--------------------------------------------------------------------------
% This program is part of the paper "...".
% Copyright (C) 2019  Alessandro De Rosis (derosis.alessandro@icloud.com)
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% This is the "D2Q9_CentralMoments_ColorGradient.m" file.
% Author: Alessandro De Rosis (derosis.alessandro@icloud.com)
% Day: 1st Feb 2020
%--------------------------------------------------------------------------

%clear all
clc
addpath('/Applications/Matlab_R2019b.app/toolbox/symbolic/symbolic')

syms U V R omega f0 f1 f2 f3 f4 f5 f6 f7 f8 Fx Fy abar...
     nu gradxR gradyR real
%% Define lattice directions, weights and other useful quantities of the D2Q9 model
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];
w = [4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.];
cs = 1/sqrt(3);
cs2 = cs^2;
cs3 = cs^3;
cs4 = cs^4;
cs6 = cs^6;
cs8 = cs^8;

f = [f0 f1 f2 f3 f4 f5 f6 f7 f8]'; % hydrodynamic populations
feq = sym(zeros(9,1));

Force = sym(zeros(9,1)); % generic force vector
T = sym(zeros(9,9)); %transformation matrix for central moments
M = zeros(9,9);  %transformation matrix for raw moments
Lambda = diag([1, 1, 1, 1, omega, omega, 1, 1, 1]);
Id = eye(9,9); % identity matrix

phi = sym(zeros(9,1));
varphi = sym(zeros(9,1));
psi = sym(zeros(9,1));
xi = sym(zeros(9,1));
correction = sym(zeros(9,1));
for i=1:length(cx)
    if(i==1)
        phi(i) = 0;
        varphi(i) = 1;
        psi(i) = -8/3;
        xi(i) = 0;
    elseif(i>1 && i<=5)
        phi(i) = 1/5;
        varphi(i) = -1/5;
        psi(i) = -1/6;
        xi(i) = 1/2;
    else
        phi(i) = 1/20;
        varphi(i) = -1/20;
        psi(i) = 1/12;
        xi(i) = 1/8;
    end
end
u = [U V];
gradR = [gradxR gradyR];
G = u'*gradR+(u'*gradR)';
for i=1:length(cx)
    % build the complete equilibria
    first_order = 1/cs2*(U*cx(i)+V*cy(i));
    second_order = 1/(2*cs4)*((cx(i)*cx(i)-1/3)*U^2+...
                              (cy(i)*cy(i)-1/3)*V^2+...
                              2*cx(i)*cy(i)*U*V);
    third_order = 1/(2*cs6)*((cx(i)^2-1/3)*cy(i)*U*U*V+(cy(i)^2-1/3)*cx(i)*U*V*V);
    fourth_order = 1/(4*cs8)*((cx(i)^2-1/3)*(cy(i)^2-1/3)*U*U*V*V);
    feq(i) = R*(phi(i)+varphi(i)*abar+w(i)*(first_order+second_order+third_order+fourth_order));
    
    ci = [cx(i) cy(i)];
    correction(i) = nu*( psi(i)*dot(u,gradR) + xi(i)*sum(sum(G*(ci'*ci))) );
    %feq(i) = feq(i) + correction(i);
    % build the complete forcing terms
    hat_cx = cx(i)/cs;
    hat_cy = cy(i)/cs;
    hat_ = [hat_cx, hat_cy];
    for a=1:2
        H1(a) = hat_(a);
        for b=1:2
            H2(a,b) = hat_(a)*hat_(b)-Id(a,b);
            for c=1:2
                hat_I = hat_(a)*Id(b,c)+hat_(b)*Id(a,c)+hat_(c)*Id(a,b);
                H3(a,b,c) = hat_(a)*hat_(b)*hat_(c)-hat_I;
                for d=1:2
                    hat_II = hat_(a)*hat_(b)*Id(c,d)+hat_(a)*hat_(c)*Id(b,d)+...
                             hat_(a)*hat_(d)*Id(b,c)+hat_(b)*hat_(c)*Id(a,d)+...
                             hat_(b)*hat_(d)*Id(a,c)+hat_(c)*hat_(d)*Id(a,b);
                    II = Id(a,b)*Id(c,d)+Id(a,c)*Id(b,d)+Id(a,d)*Id(b,c);
                    H4(a,b,c,d) = hat_(a)*hat_(b)*hat_(c)*hat_(d)-hat_II+II;
                end
            end
        end
    end
    first_order = 1/cs*(Fx*H1(1)+Fy*H1(2));                       
    second_order = 1/(2*cs2)*( (Fx*U+U*Fx)*H2(1,1) +...
                               (Fy*V+V*Fy)*H2(2,2) +...
                               (Fx*V+Fy*U)*H2(1,2) +...
                               (Fy*U+Fx*V)*H2(2,1) );
    third_order =  1/(6*cs3)*( (Fx*V*V+U*Fy*V+U*V*Fy)*H3(1,2,2) +...
                               (Fy*U*V+V*Fx*V+V*U*Fy)*H3(2,1,2) +... 
                               (Fy*V*U+V*Fy*U+V*V*Fx)*H3(2,2,1) +...
                               (Fx*U*V+U*Fx*V+U*U*Fy)*H3(1,1,2) +...
                               (Fy*U*U+V*Fx*U+V*U*Fx)*H3(2,1,1) +...
                               (Fx*V*U+U*Fy*U+U*V*Fx)*H3(1,2,1) );
    fourth_order = 1/(24*cs4)*( (Fx*U*V*V+U*Fx*V*V+U*U*Fy*V+U*U*V*Fy)*H4(1,1,2,2)+...
                                (Fx*V*U*V+U*Fy*U*V+U*V*Fx*V+U*V*U*Fy)*H4(1,2,1,2)+...
                                (Fx*V*V*U+U*Fy*V*U+U*V*Fy*U+U*V*V*Fx)*H4(1,2,2,1)+...
                                (Fy*U*U*V+V*Fx*U*V+V*U*Fx*V+V*U*U*Fy)*H4(2,1,1,2)+...
                                (Fy*U*V*U+V*Fx*V*U+V*U*Fy*U+V*U*V*Fx)*H4(2,1,2,1)+...
                                (Fy*V*U*U+V*Fy*U*U+V*V*Fx*U+V*V*U*Fx)*H4(2,2,1,1) );
    Force(i) = w(i)*(first_order + second_order + third_order + fourth_order)/R;
    
    % build the transformation matrix T 
    CX = cx(i)-U;
    CY = cy(i)-V;
    T(1,i) = 1;
    T(2,i) = CX;
    T(3,i) = CY;
    T(4,i) = CX*CX+CY*CY;
    T(5,i) = CX*CX-CY*CY;
    T(6,i) = CX*CY;
    T(7,i) = CX*CX*CY;
    T(8,i) = CX*CY*CY;
    T(9,i) = CX*CX*CY*CY;
    
    % build the tranformation matrix M
    CX = cx(i);
    CY = cy(i);
    M(1,i) = 1;
    M(2,i) = CX;
    M(3,i) = CY;
    M(4,i) = CX*CX+CY*CY;
    M(5,i) = CX*CX-CY*CY;
    M(6,i) = CX*CY;
    M(7,i) = CX*CX*CY;
    M(8,i) = CX*CY*CY;
    M(9,i) = CX*CX*CY*CY;
 end
T = simplify(T);
T = M;
%N = simplify(T*M^(-1)); %shift matrix
N = Id;
syms k0_pre k1_pre k2_pre k3_pre k4_pre k5_pre k6_pre k7_pre k8_pre real
syms k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star real
R_pre = simplify(M*f)
K_pre = simplify(N*R_pre)
R_eq = simplify(M*feq);
K_eq = simplify(N*R_eq) 
R_force = simplify(M*Force);
K_force = simplify(N*R_force); % forcing term central moments
K_pre(5) = k4_pre;
K_pre(6) = k5_pre;
K_star = simplify((Id-Lambda)*K_pre + Lambda*K_eq + (Id-Lambda/2)*K_force ) %post-collision central moments
%post-collision populations
K_sym = [R k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star];
for i=1:9
    if(K_star(i)~=sym(0))
        K_star(i) = K_sym(i);
    end
end
f_post_collision_onestep = collect(simplify(T \ K_star), K_star);
% two-steps approach
raw_moments = simplify(N^(-1)*K_star)
syms r0 r1 r2 r3 r4 r5 r6 r7 r8 real
r = [r0 r1 r2 r3 r4 r5 r6 r7 r8]'; %symbolic raw moments
f_post_collision_twosteps = collect(simplify(M\r),K_star)
