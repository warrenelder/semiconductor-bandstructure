
clc;
clear all;

syms kx ky kz;
syms eyz ezx exy exx eyy ezz; 

% e2_2  = (1/2)*(exx+2i*exy-eyy);
% e2_b2 = (1/2)*(exx-2i*exy-eyy);
% e2_1  = -(ezx+1i*eyz);
% e2_b1 = (ezx-1i*eyz);
% e2_0 = sqrt(3/2)*ezz;
% e0_0=exx+eyy+ezz;
e2_2=(exx^2+2*1i*exy-eyy^2);
e2_1=-2*(ezx+1i*eyz);
e2_0=sqrt(2/3)*(2*ezz^2-exx^2-eyy^2);
e2_b1=2*(ezx-1i*eyz);
e2_b2=(exx^2-2*1i*exy-eyy^2);
e0_0 = exx+eyy+ezz;

k1_1=-(kx+1i*ky);
k1_0=sqrt(2)*kz;
k1_b1=(kx-1i*ky);

% k_i*e_jk transforms as D0, D1, D2, D3 tensors
% J=1
disp('J = 1, M = 1, 0, -1')
U1_1  = simplify(sqrt(3/5)*e2_2*k1_b1  - sqrt(3/10)*e2_1*k1_0  + sqrt(1/10)*e2_0*k1_1);
U1_0  = simplify(sqrt(3/10)*e2_1*k1_b1 - sqrt(2/5)*e2_0*k1_0   + sqrt(3/10)*e2_b1*k1_1);
U1_b1 = simplify(sqrt(3/5)*e2_b2*k1_b1 - sqrt(3/10)*e2_b1*k1_0 + sqrt(1/10)*e2_0*k1_1);
disp(U1_1)
disp(U1_0)
disp(U1_b1)
% J=2
disp('J = 2, M = 2, 1, 0, -1, -2')
U2_2  = simplify(sqrt(2/3)*e2_2*k1_0   - sqrt(1/3)*e2_1*k1_1);
U2_1  = simplify(sqrt(1/3)*e2_2*k1_b1  + sqrt(1/6)*e2_1*k1_0  - sqrt(1/2)*e2_0*k1_1);
U2_0  = simplify(sqrt(1/2)*e2_1*k1_b1  - sqrt(1/2)*e2_b1*k1_1);
U2_b1 = simplify(-sqrt(1/3)*e2_b2*k1_1 - sqrt(1/6)*e2_b1*k1_0 + sqrt(1/2)*e2_0*k1_b1);
U2_b2 = simplify(-sqrt(2/3)*e2_b2*k1_0 + sqrt(1/3)*e2_b1*k1_0);
disp(U2_2)
disp(U2_1)
disp(U2_0)
disp(U2_b1)
disp(U2_b2)
% J=3
disp('J = 3, M = 3, 2, 1, 0, -1, -2, -3')
U3_3  = simplify(e2_2*k1_1);
U3_2  = simplify(sqrt(1/3)*e2_2*k1_0    + sqrt(2/3)*e2_1*k1_1);
U3_1  = simplify(sqrt(1/15)*e2_2*k1_b1  + sqrt(8/15)*e2_1*k1_0  + sqrt(2/5)*e2_0*k1_1);
U3_0  = simplify(sqrt(1/5)*e2_1*k1_b1   + sqrt(3/5)*e2_0*k1_0   + sqrt(1/5)*e2_b1*k1_1);
U3_b1 = simplify(-sqrt(1/15)*e2_b2*k1_1 + sqrt(8/15)*e2_b1*k1_0 + sqrt(2/5)*e2_0*k1_b1);
U3_b2 = simplify(-sqrt(1/3)*e2_b2*k1_0  + sqrt(2/3)*e2_b1*k1_b1);
U3_b3 = simplify(e2_b2*k1_b1);
disp(U3_3)
disp(U3_2)
disp(U3_1)
disp(U3_0)
disp(U3_b1)
disp(U3_b2)
disp(U3_b3)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% 
% 
% 
% 
% 
% 
% 
% %M=kron(e2_2,k1_1)
% A = simplify((sqrt(1/5)*e2_2*k1_1));
% B = simplify((-sqrt(1/10)*e2_1*k1_1));
% C = simplify((sqrt(1/30)*e2_0*k1_1));
% 
% D = simplify((-sqrt(1/10)*e2_1*k1_0));
% E = simplify((sqrt(2/15)*e2_0*k1_0));
% F = simplify((-sqrt(1/10)*e2_b1*k1_0));
% 
% G = simplify((sqrt(1/30)*e2_0*k1_b1));
% H = simplify((-sqrt(1/10)*e2_b1*k1_b1));
% I = simplify((sqrt(1/5)*e2_b2*k1_b1));
% 
% J = simplify((sqrt(1/3)*e0_0*k1_1));
% K = simplify((-sqrt(1/3)*e0_0*k1_0));
% L = simplify((sqrt(1/3)*e0_0*k1_b1));
% 
% Const=sqrt(10);
% 
% X1=simplify((A+I)*(Const/sqrt(2)));
% disp(A)
% disp(B)
% disp(C)
% disp(D)
% disp(E)
% disp(F)
% disp(G)
% disp(H)
% disp(I)
% disp(J)
% disp(K)
% disp(L)
% 
% 
% e0_0=exx+eyy+ezz;
% e2_2=(exx^2+2*1i*exy-eyy^2);
% e2_1=-2*(ezx+1i*eyz);
% e2_0=sqrt(2/3)*(2*ezz^2-exx^2-eyy^2);
% e2_b1=2*(ezx-1i*eyz);
% e2_b2=(exx^2-2*1i*exy-eyy^2);
% 
% k1_1=-(kx+1i*ky);
% k1_0=sqrt(2)*kz;
% k1_b1=(kx-1i*ky);
% 
% ek1_1=-sqrt(1/10)*e2_1*k1_0+sqrt(1/30)*e2_0*k1_1;
% ek1_0=-sqrt(2/15)*e2_0*k1_0-sqrt(1/10)*e2_1*k1_b1-sqrt(1/10)*e2_b1*k1_1;
% ek1_b1=sqrt(1/30)*e2_0*k1_b1-sqrt(1/10)*e2_b1*k1_0;
% 
% disp('From l=1, Gamma 4-:');
% ek1_const=15/sqrt(10);
% simplify((ek1_b1-ek1_1)/sqrt(2)*ek1_const)
% simplify(1i*(ek1_b1+ek1_1)/sqrt(2)*ek1_const)
% simplify(ek1_0*ek1_const)
% 
% ek2_2=-sqrt(2/15)*e2_2*k1_0+sqrt(1/15)*e2_1*k1_1;
% ek2_1=sqrt(1/15)*e2_2*k1_b1+sqrt(1/30)*e2_1*k1_0-sqrt(1/10)*e2_0*k1_1;
% ek2_0=sqrt(1/10)*e2_1*k1_b1+sqrt(1/10)*e2_b1*k1_1;
% ek2_b1=sqrt(1/10)*e2_0*k1_b1-sqrt(1/30)*e2_b1*k1_0-sqrt(1/15)*e2_b2*k1_1;
% ek2_b2=sqrt(2/15)*e2_b2*k1_0-sqrt(1/15)*e2_b1*k1_b1;
% 
% disp('From l=2, Gamma 1+: not correct');
% ek2_g1_const=5/2/sqrt(10);
% simplify((sqrt(3/2)*(ek2_b2+ek2_2)+ek2_0)*ek2_g1_const)
% 
% disp('From l=2, Gamma 3+:');
% ek2_g3_const=5/2/sqrt(10);
% simplify(sqrt(1/2)*(ek2_b2+ek2_2)*ek2_g3_const)
% simplify(ek2_0*ek2_g3_const)
% 
% disp('From l=2, Gamma_5+:');
% ek2_g5_const=sqrt(15/2);
% simplify(ek2_g5_const*1i*(ek2_b1+ek2_1)/sqrt(2))
% simplify(ek2_g5_const*(ek2_b1-ek2_1)/sqrt(2))
% simplify(ek2_g5_const*1i*(ek2_b2-ek2_2)/sqrt(2))
% 
% ek3_3=sqrt(1/7)*e2_2*k1_1;
% ek3_2=-sqrt(1/21)*e2_2*k1_0-sqrt(2/21)*e2_1*k1_1;
% ek3_1=sqrt(1/105)*e2_2*k1_b1+sqrt(8/105)*e2_1*k1_0+sqrt(2/35)*e2_0*k1_1;
% ek3_0=-sqrt(2/5)*e2_0*k1_0+sqrt(6/20)*e2_1*k1_b1+sqrt(6/20)*e2_b1*k1_1;
% ek3_b1=sqrt(2/35)*e2_0*k1_b1+sqrt(8/105)*e2_b1*k1_0+sqrt(1/105)*e2_b2*k1_1;
% ek3_b2=-sqrt(1/21)*e2_b2*k1_0-sqrt(2/21)*e2_b1*k1_b1;
% ek3_b3=sqrt(1/10)*e2_b2*k1_b1;
% 
% disp('From l=2, Gamma_5-:');
% ek3_g4p_const=1;
% simplify(ek3_g4p_const*(ek3_b3-ek3_3)/sqrt(2))
% simplify(ek3_g4p_const*1i*(ek3_b3+ek3_3)/sqrt(2))
% simplify(ek3_g4p_const*(ek3_b2+ek3_2)/sqrt(2))
