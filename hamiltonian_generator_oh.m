 clc 
 clear all
% Wigner_Eckart calculation returns interactions between near set j=1/2,
% j=3/2, j=5/2 with remote states j=1/2, j=3/2, j=5/2, j=7/2 giving a set 
% of matrices. Then using the orbital symmetrised basis functions a set of
% matrices may be derived for all possible couplings of near j states with
% remote j states

%%%% Material system codes are as follows;

%%%% <1/2,mj|1/2,mj>,           the code is M1.
%%%% <1/2,mj|3/2,mj>,           the code is M2.

%%%% <3/2,mj|1/2,mj>,           the code is M3.
%%%% <3/2,mj|3/2,mj>,           the code is M4.
%%%% <3/2,mj|5/2,mj>,           the code is M5.

%%%% <5/2,mj|3/2,mj>,           the code is M6.
%%%% <5/2,mj|5/2,mj>,           the code is M7.
%%%% <5/2,mj|7/2,mj>,           the code is M8.

Matrix = 'M7';

flag_1 = 3;
% Define orbital symmetrised Basis from; Onodera et al, Phys. Soc Jap, 21, 11
% (1966)

R_g6p_sj1_2=sym([ 1, 0; 
                  0, 1; ]);

R_g6m_pj1_2=(R_g6p_sj1_2)';

R_g7p_dj5_2=sym([ (-sqrt(1/6)), 0, 0, 0, (sqrt(5/6)), 0;
                  0, (sqrt(5/6)), 0, 0, 0, (-sqrt(1/6));]); 

R_g7m_fj5_2=(R_g7p_dj5_2)';
              
R_g8p_dj3_2=sym([ 1, 0, 0, 0; 
                  0, 1, 0, 0;
                  0, 0, 1, 0;
                  0, 0, 0, 1; ]);   
              
R_g8m_pj3_2=(R_g8p_dj3_2)';

R_g8m_fj5_2=sym([ 0, (-sqrt(1/6)),  0, 0,  0, (-sqrt(5/6));
                  0, 0,             1, 0,  0, 0;
                  0, 0,             0, -1, 0, 0;
                  (sqrt(5/6)), 0,   0, 0, (sqrt(1/6)), 0;]);

R_g6m_fj7_2=sym([ 0, 0,  0, (sqrt(7/12)),  0, 0, 0, (sqrt(5/12));
                  (-sqrt(5/12)), 0,  0, 0,  (-sqrt(7/12)), 0, 0, 0;]);
              
R_g8p_f7_2=sym([ 0,              0,              (sqrt(7/12)), 	0; 
                0,              0,              0,              (sqrt(1/4));
                (sqrt(3/4)),    0,              0,              0;
                0,              (-sqrt(5/12)),  0,              0;
                0,              0,              (-sqrt(5/12)),  0; 
                0,              0,              0,              (sqrt(3/4));
                (sqrt(1/4)),    0,              0,              0;
                0,              (sqrt(7/12)),   0,              0; ]);

% Higher order basis
R_g8_d5_2=sym([ 0,           (-sqrt(1/6)), 0, 0,  0,           (-sqrt(5/6)); 
                0,           0,            1, 0,  0,           0;
                0,           0,            0, -1, 0,           0;
                (sqrt(5/6)), 0,            0, 0,  (sqrt(1/6)), 0; ]);
          
R_g6_f7_2=sym([ 0, (-sqrt(5/12)); 
                0, 0;
                0, 0;
                (sqrt(7/12)), 0;
                0, (-sqrt(7/12)); 
                0, 0;
                0, 0;
                (sqrt(5/12)), 0; ]);

R_g7_f5_2=sym([ (-sqrt(1/6)), 0;
                0, (sqrt(5/6));
                0, 0;
                0, 0; 
                (sqrt(5/6)), 0;
                0, (-sqrt(1/6)); ]);            
            
R_g7_f7_2=sym([ 0, 0; 
                (-sqrt(3/4)), 0;
                0, (-sqrt(1/4));
                0, 0;
                0, 0; 
                (sqrt(1/4)), 0;
                0, (sqrt(3/4));
                0, 0; ]);

R_g8_f5_2=sym([ 0,              0,              0,              (sqrt(5/6));
                (-sqrt(1/6)),   0,              0,              0;
                0,              (1),            0,              0;
                0,              0,              (-1),           0; 
                0,              0,              0,              (sqrt(1/6));
                (-sqrt(5/6)),   0,              0,              0; ]);
            
R_g8_f7_2=sym([ 0,              0,              (sqrt(7/12)), 	0; 
                0,              0,              0,              (sqrt(1/4));
                (sqrt(3/4)),    0,              0,              0;
                0,              (-sqrt(5/12)),  0,              0;
                0,              0,              (-sqrt(5/12)),  0; 
                0,              0,              0,              (sqrt(3/4));
                (sqrt(1/4)),    0,              0,              0;
                0,              (sqrt(7/12)),   0,              0; ]);


switch Matrix
case 'M1'
%<1/2,mj|1/2,mj>
N_x=sym([ (0), (1/(sqrt(6)));
          (1/(sqrt(6))), (0);]);
   
N_y=sym([ (0), (-1i/(sqrt(6)));
          (1i/(sqrt(6))), (0);]);
      
N_z=sym([(1/(sqrt(6))), (0);
         (0), (-1/(sqrt(6)));]);
      
%M matrix for <Gamma_6m|Gamma_6p>
disp('<Gamma_6p|Gamma_6m>')
M_g6p_g6m_x=simplify((R_g6p_sj1_2*N_x)*R_g6m_pj1_2);
M_g6p_g6m_y=simplify((R_g6p_sj1_2*N_y)*R_g6m_pj1_2);
M_g6p_g6m_z=simplify((R_g6p_sj1_2*N_z)*R_g6m_pj1_2);

M_x=simplify((R_g6p_sj1_2*N_x)*R_g6m_pj1_2)
M_y=simplify((R_g6p_sj1_2*N_y)*R_g6m_pj1_2)
M_z=simplify((R_g6p_sj1_2*N_z)*R_g6m_pj1_2)

case 'M2'
%<1/2,mj|3/2,mj>
N_x=sym([ ((1/(2*sqrt(2)))), (0), (-(1/(2*sqrt(6)))), (0);
          (0), ((1/(2*sqrt(6)))), (0), (-(1/(2*sqrt(2))));]);
   
N_y=sym([ ((1i/(2*sqrt(2)))), (0), ((1i/(2*sqrt(6)))), (0);
          (0), ((1i/(2*sqrt(6)))), (0), (-(1i/(2*sqrt(2))));]);
      
N_z=sym([(0), (-(1/(sqrt(6)))), (0), 0;
          (0), (0), (-(1/(sqrt(6)))), (0);]);

disp('<Gamma_6p|Gamma_8m>')
M_g6p_g8m_x=simplify((R_g6p_sj1_2*-N_x)*R_g8m_pj3_2);
M_g6p_g8m_y=simplify((R_g6p_sj1_2*-N_y)*R_g8m_pj3_2);
M_g6p_g8m_z=simplify((R_g6p_sj1_2*-N_z)*R_g8m_pj3_2);

M_x=simplify((R_g6p_sj1_2*-N_x)*R_g8m_pj3_2)
M_y=simplify((R_g6p_sj1_2*-N_y)*R_g8m_pj3_2)
M_z=simplify((R_g6p_sj1_2*-N_z)*R_g8m_pj3_2)


case 'M3'
%<3/2,mj|1/2,mj>
N_x=sym([ (-(1/(2*sqrt(2)))), (0);
          (0), (-(1/(2*sqrt(6))));
          ((1/(2*sqrt(6)))), (0);
          (0), ((1/(2*sqrt(2))));]);
   
N_y=sym([ ((1i/(2*sqrt(2)))), (0);
          (0), ((1i/(2*sqrt(6))));
          ((1i/(2*sqrt(6)))), (0);
          (0), ((1i/(2*sqrt(2))));]);
      
N_z=sym([ (0), (0);
          ((1/(sqrt(6)))), (0);
          (0), ((1/(sqrt(6))));
          (0), (0);]);
      
%M matrix for <Gamma_8p|Gamma_6m>
% disp('<Gamma_8p|Gamma_6m>')
% 
% M_g8p_g6m_x=simplify((R_g8p_dj3_2*N_x)*R_g6m_pj1_2);
% M_g8p_g6m_y=simplify((R_g8p_dj3_2*N_y)*R_g6m_pj1_2);
% M_g8p_g6m_z=simplify((R_g8p_dj3_2*N_z)*R_g6m_pj1_2);
% 
% M_x=simplify((R_g8p_dj3_2*N_x)*R_g6m_pj1_2);
% M_y=simplify((R_g8p_dj3_2*N_y)*R_g6m_pj1_2);
% M_z=simplify((R_g8p_dj3_2*N_z)*R_g6m_pj1_2);

%M matrix for <Gamma_8m|Gamma_6p>
disp('<Gamma_8m|Gamma_6p>')

M_g8m_g6p_x=simplify((R_g8m_pj3_2*N_x)*R_g6p_sj1_2);
M_g8m_g6p_y=simplify((R_g8m_pj3_2*N_y)*R_g6p_sj1_2);
M_g8m_g6p_z=simplify((R_g8m_pj3_2*N_z)*R_g6p_sj1_2);

M_x=simplify((R_g8m_pj3_2*N_x)*R_g6p_sj1_2)
M_y=simplify((R_g8m_pj3_2*N_y)*R_g6p_sj1_2)
M_z=simplify((R_g8m_pj3_2*N_z)*R_g6p_sj1_2)

case 'M4'
%<3/2,mj|3/2,mj>
N_x=sym([ (0), (1/(2*sqrt(5))), (0), (0);
          (1/(2*sqrt(5))), (0), (1/(sqrt(15))), (0);
          (0), (1/(sqrt(15))), (0), (1/(2*sqrt(5)));
          (0), (0), (1/(2*sqrt(5))), (0);]);
 
N_y=sym([ (0), (-1i/(2*sqrt(5))), (0), (0);
          (1i/(2*sqrt(5))), (0), (-1i/(sqrt(15))), (0);
          (0), (1i/(sqrt(15))), (0), (-1i/(2*sqrt(5)));
          (0), (0), (1i/(2*sqrt(5))), (0);]);
     
N_z=sym([(3/(2*sqrt(15))), (0), (0), 0;
          (0), (1/(2*sqrt(15))), (0), (0);
          (0), (0), (-1/(2*sqrt(15))), (0);
          (0), (0), (0), (-3/(2*sqrt(15)));]);

%M matrix for <Gamma_8p|Gamma_8m>
disp('<Gamma_8p|Gamma_8m_3_2>')
M_g8p_g8m_x=simplify((R_g8p_dj3_2*N_x)*R_g8m_pj3_2);
M_g8p_g8m_y=simplify((R_g8p_dj3_2*N_y)*R_g8m_pj3_2);
M_g8p_g8m_z=simplify((R_g8p_dj3_2*N_z)*R_g8m_pj3_2);

M_x=simplify((R_g8p_dj3_2*N_x)*R_g8m_pj3_2);
M_y=simplify((R_g8p_dj3_2*N_y)*R_g8m_pj3_2);
M_z=simplify((R_g8p_dj3_2*N_z)*R_g8m_pj3_2);

case 'M5'
%<3/2,mj|5/2,mj>
N_x=sym([(1/(2*sqrt(3))), (0), (-1/(2*sqrt(30))), 0, 0, 0;
          (0), (1/(2*sqrt(5))), (0), (-1/(2*sqrt(10))), 0, 0;
          (0), (0), (1/(2*sqrt(10))), (0), (-1/(2*sqrt(5))), 0;
          (0), (0), (0), (1/(2*sqrt(30))), (0), (-1/(2*sqrt(3)));]);
   
N_y=sym([(1i/(2*sqrt(3))), (0), (1i/(2*sqrt(30))), 0, 0, 0;
          (0), (1i/(2*sqrt(5))), (0), (1i/(2*sqrt(10))), 0, 0;
          (0), (0), (1i/(2*sqrt(10))), (0), (1i/(2*sqrt(5))), 0;
          (0), (0), (0), (1i/(2*sqrt(30))), (0), (1i/(2*sqrt(3)));]);
      
N_z=sym([(0), ((-1/(sqrt(15)))), (0), (0), (0), (0);
          (0), (0), ((-1/(sqrt(10)))), (0), (0), (0);
          (0), (0), (0), ((-1/(sqrt(10)))), (0), (0);
          (0), (0), (0), (0), ((-1/(sqrt(15)))), (0);]);
      
disp('<Gamma_8p|Gamma_8m_5_2>')
M_g8p_g8m_x=simplify((R_g8p_dj3_2*-N_x)*R_g7p_dj5_2');
M_g8p_g8m_y=simplify((R_g8p_dj3_2*-N_y)*R_g7p_dj5_2');
M_g8p_g8m_z=simplify((R_g8p_dj3_2*-N_z)*R_g7p_dj5_2');

M_x=simplify((R_g8p_dj3_2*-N_x)*R_g7p_dj5_2')
M_y=simplify((R_g8p_dj3_2*-N_y)*R_g7p_dj5_2')
M_z=simplify((R_g8p_dj3_2*-N_z)*R_g7p_dj5_2')

% disp('<Gamma_8p|Gamma_7m_5_2>')
% M_g8p_g7m_x=simplify((R_g8p_dj3_2*-N_x)*R_g7m_fj5_2);
% M_g8p_g7m_y=simplify((R_g8p_dj3_2*-N_y)*R_g7m_fj5_2);
% M_g8p_g7m_z=simplify((R_g8p_dj3_2*-N_z)*R_g7m_fj5_2);
% 
% M_x1=simplify((R_g8p_dj3_2*-N_x)*R_g7m_fj5_2)
% M_y1=simplify((R_g8p_dj3_2*-N_y)*R_g7m_fj5_2)
% M_z1=simplify((R_g8p_dj3_2*-N_z)*R_g7m_fj5_2)   
 
case 'M6'
%<5/2,mj|3/2,mj>
N_x=sym([ (-1/(2*sqrt(3))), (0),            (0),             (0);
          (0),              (-1/(sqrt(20))), (0),            (0);
          (1/(sqrt(120))),  (0),            (-1/(sqrt(40))), (0);
          (0),              (1/(sqrt(40))), (0),             (-1/(sqrt(120)));
          (0),              (0),            (1/(sqrt(20))),  (0);
          (0),              (0),            (0),             (1/(2*sqrt(3)));]);
   
N_y=sym([ (1i/(2*sqrt(3))), (0), (0), (0);
          (0), (1i/(sqrt(20))), (0), (0);
          (1i/(sqrt(120))), (0), (1i/(sqrt(40))), (0);
          (0), (1i/(sqrt(40))), (0), (1i/(sqrt(120)));
          (0), (0), (1i/(sqrt(20))), (0);
          (0), (0), (0), (1i/(2*sqrt(3)));]);
      
N_z=sym([ (0), (0), (0), (0);
          (1/(sqrt(15))), (0), (0), (0);
          (0), (1/(sqrt(10))), (0), (0);
          (0), (0), (1/(sqrt(10))), (0);
          (0), (0), (0), (1/(sqrt(15)));
          (0), (0), (0), (0);]);
      
disp('<Gamma_7m|Gamma_8p_3_2>')
M_g7m_g8p_x=simplify((R_g7p_dj5_2*N_x)*R_g8m_pj3_2);
M_g7m_g8p_y=simplify((R_g7p_dj5_2*N_y)*R_g8m_pj3_2);
M_g7m_g8p_z=simplify((R_g7p_dj5_2*N_z)*R_g8m_pj3_2);

M_x=simplify((R_g7p_dj5_2*N_x)*R_g8m_pj3_2);
M_y=simplify((R_g7p_dj5_2*N_y)*R_g8m_pj3_2);
M_z=simplify((R_g7p_dj5_2*N_z)*R_g8m_pj3_2);      
      
% disp('<Gamma_8m_5_2|Gamma_8p>')
% M_g8p_g8m_x=simplify((R_g8p_dj3_2*N_x)*R_g8m_fj5_2);
% M_g8p_g8m_y=simplify((R_g8p_dj3_2*N_y)*R_g8m_fj5_2);
% M_g8p_g8m_z=simplify((R_g8p_dj3_2*N_z)*R_g8m_fj5_2);

% M_x=simplify((R_g8m_fj5_2*N_x)*R_g8p_dj3_2)
% M_y=simplify((R_g8m_fj5_2*N_y)*R_g8p_dj3_2)
% M_z=simplify((R_g8m_fj5_2*N_z)*R_g8p_dj3_2)

 disp('<Gamma_7m|Gamma_8m_3_2>')
% M_g7m_g8p_x=simplify((R_g7m_fj5_2'*N_x)*R_g8p_dj3_2);
% M_g7m_g8p_y=simplify((R_g7m_fj5_2'*N_y)*R_g8p_dj3_2);
% M_g7m_g8p_z=simplify((R_g7m_fj5_2'*N_z)*R_g8p_dj3_2);
% 
M_x=simplify((R_g7m_fj5_2'*N_x)*R_g8p_dj3_2);
M_y=simplify((R_g7m_fj5_2'*N_y)*R_g8p_dj3_2);
M_z=simplify((R_g7m_fj5_2'*N_z)*R_g8p_dj3_2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'M7'
%<5/2,mj|5/2,mj>
N_x=sym([ (0), (-sqrt(2/42)*-sqrt(1/2)), 0, 0, 0, 0;
          (sqrt(2/42)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0, 0, 0;
          (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(18/210)*-sqrt(1/2)), 0, 0;
          (0), (0), (sqrt(18/210)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0;
          (0), (0), (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(2/42)*-sqrt(1/2));
          (0), (0), (0), (0),(sqrt(2/42)*sqrt(1/2)), (0);]);
   
N_y=sym([ (0), (-1i*sqrt(2/42)*sqrt(1/2)), 0, 0, 0, 0;
          (1i*sqrt(2/42)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0, 0, 0;
          (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(18/210)*sqrt(1/2)), 0, 0;
          (0), (0), (1i*sqrt(18/210)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0;
          (0), (0), (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(2/42)*sqrt(1/2));
          (0), (0), (0), (0),(1i*sqrt(2/42)*sqrt(1/2)), (0);]);
      
N_z=sym([ (sqrt(5/42)), (0), 0, 0, 0, 0;
          (0), (sqrt(9/210)), (0), 0, 0, 0;
          (0), (0), (sqrt(1/210)), (0), 0, 0;
          (0), (0), (0), (-sqrt(1/210)), (0), 0;
          (0), (0), (0), (0), (-sqrt(9/210)), (0);
          (0), (0), (0), (0),(0), (-sqrt(5/42));]);
      
% disp('<Gamma_7p|Gamma_7m>')
% M_g7p_g7m_x=simplify((R_g7p_dj5_2*N_x)*R_g7m_fj5_2);
% M_g7p_g7m_y=simplify((R_g7p_dj5_2*N_y)*R_g7m_fj5_2);
% M_g7p_g7m_z=simplify((R_g7p_dj5_2*N_z)*R_g7m_fj5_2);
% 
% M_x=simplify((R_g7p_dj5_2*N_x)*R_g7m_fj5_2);
% M_y=simplify((R_g7p_dj5_2*N_y)*R_g7m_fj5_2);
% M_z=simplify((R_g7p_dj5_2*N_z)*R_g7m_fj5_2);

disp('<Gamma_7p|Gamma_8m_5_2>')

M_x=simplify((R_g7p_dj5_2*N_x)*R_g8m_fj5_2');
M_y=simplify((R_g7p_dj5_2*N_y)*R_g8m_fj5_2');
M_z=simplify((R_g7p_dj5_2*N_z)*R_g8m_fj5_2');

      
case 'M8'
%<5/2,mj|7/2,mj>
N_x=sym([(sqrt(3/24)*sqrt(1/2)), (0), (sqrt(1/168)*-sqrt(1/2)), 0, 0, 0, 0, 0;
         (0), (sqrt(15/168)*sqrt(1/2)), (0), (sqrt(3/168)*-sqrt(1/2)), 0, 0, 0, 0;
         (0), (0), (sqrt(5/84)*sqrt(1/2)), (0), (sqrt(3/84)*-sqrt(1/2)), 0, 0, 0;
         (0), (0), (0), (sqrt(3/84)*sqrt(1/2)), (0), (sqrt(5/84)*-sqrt(1/2)), 0, 0;
         (0), (0), (0), (0), (sqrt(3/168)*sqrt(1/2)), (0), (sqrt(15/168)*-sqrt(1/2)), 0;
         (0), (0), (0), (0), (0),(sqrt(1/168)*sqrt(1/2)), (0), (sqrt(3/24)*-sqrt(1/2));]);
   
N_y=sym([(1i*sqrt(3/24)*sqrt(1/2)), (0), (1i*sqrt(1/168)*sqrt(1/2)), 0, 0, 0, 0, 0;
       (0), (1i*sqrt(15/168)*sqrt(1/2)), (0), (1i*sqrt(3/168)*sqrt(1/2)), 0, 0, 0, 0;
       (0), (0), (1i*sqrt(5/84)*sqrt(1/2)), (0), (1i*sqrt(3/84)*sqrt(1/2)), 0, 0, 0;
       (0), (0), (0), (1i*sqrt(3/84)*sqrt(1/2)), (0), (1i*sqrt(5/84)*sqrt(1/2)), 0, 0;
       (0), (0), (0), (0), (1i*sqrt(3/168)*sqrt(1/2)), (0), (1i*sqrt(15/168)*sqrt(1/2)), 0;
       (0), (0), (0), (0), (0),(1i*sqrt(1/168)*sqrt(1/2)), (0), (1i*sqrt(3/24)*sqrt(1/2));]);

N_z=sym([(0), (-sqrt(3/84)), (0), 0, 0, 0, 0, 0;
       (0), (0), (-sqrt(5/84)), (0), 0, 0, 0, 0;
       (0), (0), (0), (-sqrt(3/42)), (0), 0, 0, 0;
       (0), (0), (0), (0), (-sqrt(3/42)), (0), 0, 0;
       (0), (0), (0), (0), (0), (-sqrt(5/84)), (0), 0;
       (0), (0), (0), (0), (0),(0), (-sqrt(3/84)), (0);]);
   
% M matrix for <Gamma_8d_(5/2)|Gamma_6f_(7/2)>
M_x=simplify(((R_g7p_dj5_2*N_x)*R_g6_f7_2))
M_y=simplify(((R_g7p_dj5_2*N_y)*R_g6_f7_2))
M_z=simplify(((R_g7p_dj5_2*N_z)*R_g6_f7_2))

% % M matrix for <Gamma_8d_(5/2)|Gamma_7f_(7/2)>
% M_x=simplify(((R_g8_d5_2*N_x)*R_g7_f7_2));
% M_y=simplify(((R_g8_d5_2*N_y)*R_g7_f7_2));
% M_z=simplify(((R_g8_d5_2*N_z)*R_g7_f7_2));
 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 syms kx_l kx_r ky_l ky_r kx ky kz k_x k_y k_z;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % <G8p| - |G7p>
% 
% %<3/2,mj|3/2,mj>
% N_x=sym([ (0), (1/(2*sqrt(5))), (0), (0);
%           (1/(2*sqrt(5))), (0), (1/(sqrt(15))), (0);
%           (0), (1/(sqrt(15))), (0), (1/(2*sqrt(5)));
%           (0), (0), (1/(2*sqrt(5))), (0);]);
%  
% N_y=sym([ (0), (-1i/(2*sqrt(5))), (0), (0);
%           (1i/(2*sqrt(5))), (0), (-1i/(sqrt(15))), (0);
%           (0), (1i/(sqrt(15))), (0), (-1i/(2*sqrt(5)));
%           (0), (0), (1i/(2*sqrt(5))), (0);]);
%      
% N_z=sym([(3/(2*sqrt(15))), (0), (0), 0;
%           (0), (1/(2*sqrt(15))), (0), (0);
%           (0), (0), (-1/(2*sqrt(15))), (0);
%           (0), (0), (0), (-3/(2*sqrt(15)));]); 
% %<3/2,mj|5/2,mj>
% N_x1=sym([(1/(2*sqrt(3))), (0), (-1/(2*sqrt(30))), 0, 0, 0;
%           (0), (1/(2*sqrt(5))), (0), (-1/(2*sqrt(10))), 0, 0;
%           (0), (0), (1/(2*sqrt(10))), (0), (-1/(2*sqrt(5))), 0;
%           (0), (0), (0), (1/(2*sqrt(30))), (0), (-1/(2*sqrt(3)));]);
%    
% N_y1=sym([(1i/(2*sqrt(3))), (0), (1i/(2*sqrt(30))), 0, 0, 0;
%           (0), (1i/(2*sqrt(5))), (0), (1i/(2*sqrt(10))), 0, 0;
%           (0), (0), (1i/(2*sqrt(10))), (0), (1i/(2*sqrt(5))), 0;
%           (0), (0), (0), (1i/(2*sqrt(30))), (0), (1i/(2*sqrt(3)));]);
%       
% N_z1=sym([(0), ((-1/(sqrt(15)))), (0), (0), (0), (0);
%           (0), (0), ((-1/(sqrt(10)))), (0), (0), (0);
%           (0), (0), (0), ((-1/(sqrt(10)))), (0), (0);
%           (0), (0), (0), (0), ((-1/(sqrt(15)))), (0);]);
%       
% %<5/2,mj|5/2,mj>
% N_x2=sym([ (0), (-sqrt(2/42)*-sqrt(1/2)), 0, 0, 0, 0;
%           (sqrt(2/42)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0, 0, 0;
%           (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(18/210)*-sqrt(1/2)), 0, 0;
%           (0), (0), (sqrt(18/210)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0;
%           (0), (0), (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(2/42)*-sqrt(1/2));
%           (0), (0), (0), (0),(sqrt(2/42)*sqrt(1/2)), (0);]);
%    
% N_y2=sym([ (0), (-1i*sqrt(2/42)*sqrt(1/2)), 0, 0, 0, 0;
%           (1i*sqrt(2/42)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0, 0, 0;
%           (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(18/210)*sqrt(1/2)), 0, 0;
%           (0), (0), (1i*sqrt(18/210)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0;
%           (0), (0), (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(2/42)*sqrt(1/2));
%           (0), (0), (0), (0),(1i*sqrt(2/42)*sqrt(1/2)), (0);]);
%       
% N_z2=sym([ (sqrt(5/42)), (0), 0, 0, 0, 0;
%           (0), (sqrt(9/210)), (0), 0, 0, 0;
%           (0), (0), (sqrt(1/210)), (0), 0, 0;
%           (0), (0), (0), (-sqrt(1/210)), (0), 0;
%           (0), (0), (0), (0), (-sqrt(9/210)), (0);
%           (0), (0), (0), (0),(0), (-sqrt(5/42));]);
% if flag_1 == 1
% disp('<G8p| - |G7m>-<G7m| - |G7p>') 
% M_x1=simplify((R_g8p_dj3_2*-N_x1)*R_g7m_fj5_2);
% M_y1=simplify((R_g8p_dj3_2*-N_y1)*R_g7m_fj5_2);
% M_z1=simplify((R_g8p_dj3_2*-N_z1)*R_g7m_fj5_2);
%  
% M_x2=simplify((R_g7m_fj5_2'*N_x2)*R_g7m_fj5_2);
% M_y2=simplify((R_g7m_fj5_2'*N_y2)*R_g7m_fj5_2);
% M_z2=simplify((R_g7m_fj5_2'*N_z2)*R_g7m_fj5_2); 
% elseif flag_1 == 2
% disp('<G8p| - |G8m1>-<G8m1| - |G7p>') 
% M_x1=simplify((R_g8p_dj3_2*N_x)*R_g8p_dj3_2);
% M_y1=simplify((R_g8p_dj3_2*N_y)*R_g8p_dj3_2);
% M_z1=simplify((R_g8p_dj3_2*N_z)*R_g8p_dj3_2);
%  
% M_x2=simplify((R_g8p_dj3_2'*-N_x1)*R_g7m_fj5_2);
% M_y2=simplify((R_g8p_dj3_2'*-N_y1)*R_g7m_fj5_2);
% M_z2=simplify((R_g8p_dj3_2'*-N_z1)*R_g7m_fj5_2); 
% elseif flag_1 == 3
% disp('<G8p| - |G8m2>-<G8m2| - |G7p>') 
% M_x1=simplify((R_g8p_dj3_2*(-N_x1))*R_g8m_fj5_2');
% M_y1=simplify((R_g8p_dj3_2*(-N_y1))*R_g8m_fj5_2');
% M_z1=simplify((R_g8p_dj3_2*(-N_z1))*R_g8m_fj5_2');
%  
% M_x2=simplify((R_g8m_fj5_2*N_x2)*R_g7m_fj5_2);
% M_y2=simplify((R_g8m_fj5_2*N_y2)*R_g7m_fj5_2);
% M_z2=simplify((R_g8m_fj5_2*N_z2)*R_g7m_fj5_2);     
%     
% end    
% disp(M_x1) 
% disp(M_y1) 
% disp(M_z1) 
% disp(M_x2) 
% disp(M_y2) 
% disp(M_z2)
% 
% a12=simplify((k_x.*M_x1)*(M_y2.*k_y));
% b12=simplify((k_y.*M_y1)*(M_x2.*k_x));
% c12=simplify(a12+b12);
% d12=simplify((k_x.*M_x1)*(M_x2.*k_x));
% e12=simplify((k_y.*M_y1)*(M_y2.*k_y));
% f12=simplify((k_z.*M_z1)*(M_z2.*k_z));
% p12=simplify((k_x.*M_x1)*(M_z2.*k_z));
% r12=simplify((k_z.*M_z1)*(M_x2.*k_x));
% s12=simplify((k_z.*M_z1)*(M_y2.*k_y));
% q12=simplify((k_y.*M_y1)*(M_z2.*k_z));
% 
% 
% M1221=simplify(((945)/(sqrt((105))))*(c12+d12+e12+f12+p12+r12+s12+q12));
% disp('Un-Ordered Lowdin interaction matrices')
% disp(M1221)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
a=simplify((kx_l.*M_x)*(M_y'.*ky_r));
b=simplify((ky_l.*M_y)*(M_x'.*kx_r));
c=simplify(a+b);
d=simplify((kx_l.*M_x)*(M_x'.*kx_r));
e=simplify((ky_l.*M_y)*(M_y'.*ky_r));
f=simplify((k_z.*M_z)*(M_z'.*k_z));
p=simplify((kx_l.*M_x)*(M_z'.*k_z));
r=simplify((k_z.*M_z)*(M_x'.*kx_r));
s=simplify((k_z.*M_z)*(M_y'.*ky_r));
q=simplify((ky_l.*M_y)*(M_z'.*k_z));
M=c+d+e+f+p+r+s+q;
disp('Ordered Lowdin interaction matrices')
disp(M)


a1=simplify((kx.*M_x)*(M_y'.*ky));
b1=simplify((ky.*M_y)*(M_x'.*kx));
c1=simplify(a1+b1);
d1=simplify((kx.*M_x)*(M_x'.*kx));
e1=simplify((ky.*M_y)*(M_y'.*ky));
f1=simplify((kz.*M_z)*(M_z'.*kz));
p1=simplify((kx.*M_x)*(M_z'.*kz));
r1=simplify((kz.*M_z)*(M_x'.*kx));
s1=simplify((kz.*M_z)*(M_y'.*ky));
q1=simplify((ky.*M_y)*(M_z'.*kz));
M1=c1+d1+e1+f1+p1+r1+s1+q1;
disp('Un-Ordered Lowdin interaction matrices')
disp(M1)
v=simplify((kx.*M_x));
w=simplify((ky.*M_y));
z=simplify((kz.*M_z));
D1=v+w+z;
disp('Direct interaction matrices')
% disp(v)
% disp(w)
% disp(z)
 disp(D1)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a11=simplify((kx.*M_x1)*(M_y1'.*ky));
% b11=simplify((ky.*M_y1)*(M_x1'.*kx));
% c11=simplify(a11+b11);
% d11=simplify((kx.*M_x1)*(M_x1'.*kx));
% e11=simplify((ky.*M_y1)*(M_y1'.*ky));
% f11=simplify((kz.*M_z1)*(M_z1'.*kz));
% p11=simplify((kx.*M_x1)*(M_z1'.*kz));
% r11=simplify((kz.*M_z1)*(M_x1'.*kx));
% s11=simplify((kz.*M_z1)*(M_y1'.*ky));
% q11=simplify((ky.*M_y1)*(M_z1'.*kz));
% M11=c11+d11+e11+f11+p11+r11+s11+q11;
% disp('Un-Ordered Lowdin interaction matrices G7')
% disp(M11)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %<3/2,mj|3/2,mj>
% N_x1=sym([ (0), (1/(2*sqrt(5))), (0), (0);
%           (1/(2*sqrt(5))), (0), (1/(sqrt(15))), (0);
%           (0), (1/(sqrt(15))), (0), (1/(2*sqrt(5)));
%           (0), (0), (1/(2*sqrt(5))), (0);]);
%  
% N_y1=sym([ (0), (-1i/(2*sqrt(5))), (0), (0);
%           (1i/(2*sqrt(5))), (0), (-1i/(sqrt(15))), (0);
%           (0), (1i/(sqrt(15))), (0), (-1i/(2*sqrt(5)));
%           (0), (0), (1i/(2*sqrt(5))), (0);]);
%      
% N_z1=sym([(3/(2*sqrt(15))), (0), (0), 0;
%           (0), (1/(2*sqrt(15))), (0), (0);
%           (0), (0), (-1/(2*sqrt(15))), (0);
%           (0), (0), (0), (-3/(2*sqrt(15)));]);
% 
% %M matrix for <Gamma_8p|Gamma_8m>
% disp('<Gamma_8p|Gamma_8m_3_2>')
% 
% M_x1=simplify((R_g8p_dj3_2*N_x1)*R_g8m_pj3_2);
% M_y1=simplify((R_g8p_dj3_2*N_y1)*R_g8m_pj3_2);
% M_z1=simplify((R_g8p_dj3_2*N_z1)*R_g8m_pj3_2);
% disp(M_x1)
% disp(M_y1)
% disp(M_z1)
% 
% %<3/2,mj|5/2,mj>
% N_x2=sym([(1/(2*sqrt(3))), (0), (-1/(2*sqrt(30))), 0, 0, 0;
%           (0), (1/(2*sqrt(5))), (0), (-1/(2*sqrt(10))), 0, 0;
%           (0), (0), (1/(2*sqrt(10))), (0), (-1/(2*sqrt(5))), 0;
%           (0), (0), (0), (1/(2*sqrt(30))), (0), (-1/(2*sqrt(3)));]);
%    
% N_y2=sym([(1i/(2*sqrt(3))), (0), (1i/(2*sqrt(30))), 0, 0, 0;
%           (0), (1i/(2*sqrt(5))), (0), (1i/(2*sqrt(10))), 0, 0;
%           (0), (0), (1i/(2*sqrt(10))), (0), (1i/(2*sqrt(5))), 0;
%           (0), (0), (0), (1i/(2*sqrt(30))), (0), (1i/(2*sqrt(3)));]);
%       
% N_z2=sym([(0), ((-1/(sqrt(15)))), (0), (0), (0), (0);
%           (0), (0), ((-1/(sqrt(10)))), (0), (0), (0);
%           (0), (0), (0), ((-1/(sqrt(10)))), (0), (0);
%           (0), (0), (0), (0), ((-1/(sqrt(15)))), (0);]);
%       
% disp('<Gamma_8p|Gamma_8m_5_2>')
% 
% M_x2=simplify((R_g8p_dj3_2*-N_x2)*R_g8m_fj5_2');
% M_y2=simplify((R_g8p_dj3_2*-N_y2)*R_g8m_fj5_2');
% M_z2=simplify((R_g8p_dj3_2*-N_z2)*R_g8m_fj5_2');
% disp(M_x2)
% disp(M_y2)
% disp(M_z2)
% 
% a12=simplify((k_x.*M_x1)*(M_y2'.*k_y));
% b12=simplify((k_y.*M_y1)*(M_x2'.*k_x));
% c12=simplify(a12+b12);
% d12=simplify((k_x.*M_x1)*(M_x2'.*k_x));
% e12=simplify((k_y.*M_y1)*(M_y2'.*k_y));
% f12=simplify((k_z.*M_z1)*(M_z2'.*k_z));
% p12=simplify((k_x.*M_x1)*(M_z2'.*k_z));
% r12=simplify((k_z.*M_z1)*(M_x2'.*k_x));
% s12=simplify((k_z.*M_z1)*(M_y2'.*k_y));
% q12=simplify((k_y.*M_y1)*(M_z2'.*k_z));
% 
% a21=simplify((k_x.*M_x2)*(M_y1'.*k_y));
% b21=simplify((k_y.*M_y2)*(M_x1'.*k_x));
% c21=simplify(a21+b21);
% d21=simplify((k_x.*M_x2)*(M_x1'.*k_x));
% e21=simplify((k_y.*M_y2)*(M_y1'.*k_y));
% f21=simplify((k_z.*M_z2)*(M_z1'.*k_z));
% p21=simplify((k_x.*M_x2)*(M_z1'.*k_z));
% r21=simplify((k_z.*M_z2)*(M_x1'.*k_x));
% s21=simplify((k_z.*M_z2)*(M_y1'.*k_y));
% q21=simplify((k_y.*M_y2)*(M_z1'.*k_z));
% M1221=c12+d12+e12+f12+p12+r12+s12+q12+c21+d21+e21+f21+p21+r21+s21+q21;
% disp('Un-Ordered Lowdin interaction matrices')
% disp(M1221)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOWDIN interaction between VB and SO bands

%<3/2,mj|3/2,mj>
% N_x1=sym([ (0), (1/(2*sqrt(5))), (0), (0);
%           (1/(2*sqrt(5))), (0), (1/(sqrt(15))), (0);
%           (0), (1/(sqrt(15))), (0), (1/(2*sqrt(5)));
%           (0), (0), (1/(2*sqrt(5))), (0);]);
%  
% N_y1=sym([ (0), (-1i/(2*sqrt(5))), (0), (0);
%           (1i/(2*sqrt(5))), (0), (-1i/(sqrt(15))), (0);
%           (0), (1i/(sqrt(15))), (0), (-1i/(2*sqrt(5)));
%           (0), (0), (1i/(2*sqrt(5))), (0);]);
%      
% N_z1=sym([(3/(2*sqrt(15))), (0), (0), 0;
%           (0), (1/(2*sqrt(15))), (0), (0);
%           (0), (0), (-1/(2*sqrt(15))), (0);
%           (0), (0), (0), (-3/(2*sqrt(15)));]);
% 
% M matrix for <Gamma_8p|Gamma_8m>
% disp('<Gamma_8p|Gamma_8m_3_2>')
% 
% M_x1=simplify((R_g8p_dj3_2*N_x1)*R_g8m_pj3_2);
% M_y1=simplify((R_g8p_dj3_2*N_y1)*R_g8m_pj3_2);
% M_z1=simplify((R_g8p_dj3_2*N_z1)*R_g8m_pj3_2);
% disp(M_x1)
% disp(M_y1)
% disp(M_z1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% <3/2,mj|5/2,mj>
% N_x2=sym([(1/(2*sqrt(3))), (0), (-1/(2*sqrt(30))), 0, 0, 0;
%           (0), (1/(2*sqrt(5))), (0), (-1/(2*sqrt(10))), 0, 0;
%           (0), (0), (1/(2*sqrt(10))), (0), (-1/(2*sqrt(5))), 0;
%           (0), (0), (0), (1/(2*sqrt(30))), (0), (-1/(2*sqrt(3)));]);
%    
% N_y2=sym([(1i/(2*sqrt(3))), (0), (1i/(2*sqrt(30))), 0, 0, 0;
%           (0), (1i/(2*sqrt(5))), (0), (1i/(2*sqrt(10))), 0, 0;
%           (0), (0), (1i/(2*sqrt(10))), (0), (1i/(2*sqrt(5))), 0;
%           (0), (0), (0), (1i/(2*sqrt(30))), (0), (1i/(2*sqrt(3)));]);
%       
% N_z2=sym([(0), ((-1/(sqrt(15)))), (0), (0), (0), (0);
%           (0), (0), ((-1/(sqrt(10)))), (0), (0), (0);
%           (0), (0), (0), ((-1/(sqrt(10)))), (0), (0);
%           (0), (0), (0), (0), ((-1/(sqrt(15)))), (0);]);
%       
% disp('<Gamma_8p|Gamma_7m_5_2>')
% 
% M_x2=simplify((R_g8p_dj3_2*(-N_x2))*R_g7p_dj5_2')
% M_y2=simplify((R_g8p_dj3_2*(-N_y2))*R_g7p_dj5_2')
% M_z2=simplify((R_g8p_dj3_2*(-N_z2))*R_g7p_dj5_2')  %R_g8m_fj5_2'   R_g7p_dj5_2
% disp(M_x2)
% disp(M_y2)
% disp(M_z2)
% 
% disp('<Gamma_8p|Gamma_8m_5_2>')
% 
% M_x4=simplify((R_g8p_dj3_2*(-N_x2))*R_g8m_fj5_2');
% M_y4=simplify((R_g8p_dj3_2*(-N_y2))*R_g8m_fj5_2');
% M_z4=simplify((R_g8p_dj3_2*(-N_z2))*R_g8m_fj5_2');  %R_g8m_fj5_2'   R_g7p_dj5_2
% disp(M_x4)
% disp(M_y4)
% disp(M_z4)
% <5/2,mj|5/2,mj>
% N_x3=sym([ (0), (-sqrt(2/42)*-sqrt(1/2)), 0, 0, 0, 0;
%           (sqrt(2/42)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0, 0, 0;
%           (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(18/210)*-sqrt(1/2)), 0, 0;
%           (0), (0), (sqrt(18/210)*sqrt(1/2)), (0), (-sqrt(16/210)*-sqrt(1/2)), 0;
%           (0), (0), (0), (sqrt(16/210)*sqrt(1/2)), (0), (-sqrt(2/42)*-sqrt(1/2));
%           (0), (0), (0), (0),(sqrt(2/42)*sqrt(1/2)), (0);]);
%    
% N_y3=sym([ (0), (-1i*sqrt(2/42)*sqrt(1/2)), 0, 0, 0, 0;
%           (1i*sqrt(2/42)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0, 0, 0;
%           (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(18/210)*sqrt(1/2)), 0, 0;
%           (0), (0), (1i*sqrt(18/210)*sqrt(1/2)), (0), (-1i*sqrt(16/210)*sqrt(1/2)), 0;
%           (0), (0), (0), (1i*sqrt(16/210)*sqrt(1/2)), (0), (-1i*sqrt(2/42)*sqrt(1/2));
%           (0), (0), (0), (0),(1i*sqrt(2/42)*sqrt(1/2)), (0);]);
%       
% N_z3=sym([ (sqrt(5/42)), (0), 0, 0, 0, 0;
%           (0), (sqrt(9/210)), (0), 0, 0, 0;
%           (0), (0), (sqrt(1/210)), (0), 0, 0;
%           (0), (0), (0), (-sqrt(1/210)), (0), 0;
%           (0), (0), (0), (0), (-sqrt(9/210)), (0);
%           (0), (0), (0), (0),(0), (-sqrt(5/42));]);
%       
% disp('<Gamma_7m|Gamma_7p>')
% 
% M_x3=simplify((R_g7m_fj5_2'*N_x3)*R_g7p_dj5_2');
% M_y3=simplify((R_g7m_fj5_2'*N_y3)*R_g7p_dj5_2');
% M_z3=simplify((R_g7m_fj5_2'*N_z3)*R_g7p_dj5_2');
% 
% disp(M_x3)
% disp(M_y3)
% disp(M_z3)
% 
% disp('<Gamma_8m_5_2|Gamma_7p>')
% 
% M_x5=simplify((R_g8m_fj5_2*N_x3)*R_g7p_dj5_2');
% M_y5=simplify((R_g8m_fj5_2*N_y3)*R_g7p_dj5_2');
% M_z5=simplify((R_g8m_fj5_2*N_z3)*R_g7p_dj5_2');
% 
% disp(M_x5)
% disp(M_y5)
% disp(M_z5)
% 
% a12=simplify((kx.*M_x1)*(M_y2.*ky))
% b12=simplify((ky.*M_y1)*(M_x2.*kx))
% c12=simplify(a12+b12)
% d12=simplify((kx.*M_x1)*(M_x2.*kx));
% e12=simplify((ky.*M_y1)*(M_y2.*ky));
% f12=simplify((kz.*M_z1)*(M_z2.*kz));
% p12=simplify((kx.*M_x1)*(M_z2.*kz));
% r12=simplify((kz.*M_z1)*(M_x2.*kx));
% s12=simplify((kz.*M_z1)*(M_y2.*ky));
% q12=simplify((ky.*M_y1)*(M_z2.*kz));
% 
% 
% M12=c12+d12+e12+f12+p12+r12+s12+q12;
% disp('Un-Ordered Lowdin interaction matrices <Gamma_8p|Gamma_8m_3_2><Gamma_8m_3_2|Gamma_7m_5_2>')
% disp(M12)
% 
% a45=simplify((kx.*M_x4)*(M_y5.*ky))
% b45=simplify((ky.*M_y4)*(M_x5.*kx))
% c45=simplify(a45+b45)
% d45=simplify((kx.*M_x4)*(M_x5.*kx));
% e45=simplify((ky.*M_y4)*(M_y5.*ky));
% f45=simplify((kz.*M_z4)*(M_z5.*kz));
% p45=simplify((kx.*M_x4)*(M_z5.*kz));
% r45=simplify((kz.*M_z4)*(M_x5.*kx));
% s45=simplify((kz.*M_z4)*(M_y5.*ky));
% q45=simplify((ky.*M_y4)*(M_z5.*kz));
% 
% 
% M45=c45+d45+e45+f45+p45+r45+s45+q45;
% disp('Un-Ordered Lowdin interaction matrices <Gamma_8p|Gamma_8m_5_2><Gamma_8m_5_2|Gamma_7p>')
% disp(M45)
% 
% a23=simplify((kx.*M_x2)*(M_y3'.*ky))
% b23=simplify((ky.*M_y2)*(M_x3'.*kx))
% c23=simplify(a23+b23)
% d23=simplify((kx.*M_x2)*(M_x3'.*kx));
% e23=simplify((ky.*M_y2)*(M_y3'.*ky));
% f23=simplify((kz.*M_z2)*(M_z3'.*kz));
% p23=simplify((kx.*M_x2)*(M_z3'.*kz));
% r23=simplify((kz.*M_z2)*(M_x3'.*kx));
% s23=simplify((kz.*M_z2)*(M_y3'.*ky));
% q23=simplify((ky.*M_y2)*(M_z3'.*kz));
% 
% 
% M23=c23+d23+e23+f23+p23+r23+s23+q23;
% disp('Un-Ordered Lowdin interaction matrices <Gamma_8p|Gamma_7m_5_2><Gamma_7m|Gamma_7p>')
% disp(M23)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% LOWDIN interaction between VB and VV g6 bands

% %<3/2,mj|1/2,mj>
% N_x1=sym([ (-(1/(2*sqrt(2)))), (0);
%           (0), (-(1/(2*sqrt(6))));
%           ((1/(2*sqrt(6)))), (0);
%           (0), ((1/(2*sqrt(2))));]);
%    
% N_y1=sym([ ((1i/(2*sqrt(2)))), (0);
%           (0), ((1i/(2*sqrt(6))));
%           ((1i/(2*sqrt(6)))), (0);
%           (0), ((1i/(2*sqrt(2))));]);
%       
% N_z1=sym([ (0), (0);
%           ((1/(sqrt(6)))), (0);
%           (0), ((1/(sqrt(6))));
%           (0), (0);]);
% 
% %M matrix for <Gamma_8p|Gamma_8m>
% disp('<Gamma_8p|Gamma_6m>')
% 
% M_x1=simplify((R_g8p_dj3_2*N_x1)*R_g6m_pj1_2);  %R_g6m_pj1_2
% M_y1=simplify((R_g8p_dj3_2*N_y1)*R_g6m_pj1_2);
% M_z1=simplify((R_g8p_dj3_2*N_z1)*R_g6m_pj1_2);
% disp(M_x1)
% disp(M_y1)
% disp(M_z1)
% 
% %<1/2,mj|1/2,mj>
% N_x2=sym([ (0), (1/(sqrt(6)));
%           (1/(sqrt(6))), (0);]);
%    
% N_y2=sym([ (0), (-1i/(sqrt(6)));
%           (1i/(sqrt(6))), (0);]);
%       
% N_z2=sym([(1/(sqrt(6))), (0);
%          (0), (-1/(sqrt(6)));]);
%       
% disp('<Gamma_6m|Gamma_6p>.................................')
% M_x2=simplify((R_g6m_pj1_2*N_x2)*R_g6p_sj1_2);
% M_y2=simplify((R_g6m_pj1_2*N_y2)*R_g6p_sj1_2);
% M_z2=simplify((R_g6m_pj1_2*N_z2)*R_g6p_sj1_2);  
% disp(M_x2)
% disp(M_y2)
% disp(M_z2)
% 
% %<3/2,mj|3/2,mj>
% N_x3=sym([ (0), (1/(2*sqrt(5))), (0), (0);
%           (1/(2*sqrt(5))), (0), (1/(sqrt(15))), (0);
%           (0), (1/(sqrt(15))), (0), (1/(2*sqrt(5)));
%           (0), (0), (1/(2*sqrt(5))), (0);]);
%  
% N_y3=sym([ (0), (-1i/(2*sqrt(5))), (0), (0);
%           (1i/(2*sqrt(5))), (0), (-1i/(sqrt(15))), (0);
%           (0), (1i/(sqrt(15))), (0), (-1i/(2*sqrt(5)));
%           (0), (0), (1i/(2*sqrt(5))), (0);]);
%      
% N_z3=sym([(3/(2*sqrt(15))), (0), (0), 0;
%           (0), (1/(2*sqrt(15))), (0), (0);
%           (0), (0), (-1/(2*sqrt(15))), (0);
%           (0), (0), (0), (-3/(2*sqrt(15)));]);
% 
% disp('<Gamma_8p|Gamma_8m>')
% M_x3=simplify((R_g8p_dj3_2*N_x3)*R_g8m_pj3_2);
% M_y3=simplify((R_g8p_dj3_2*N_y3)*R_g8m_pj3_2);
% M_z3=simplify((R_g8p_dj3_2*N_z3)*R_g8m_pj3_2);  
% disp(M_x3)
% disp(M_y3)
% disp(M_z3)
% a12=simplify((kx.*M_x1)*(M_y2.*ky));
% b12=simplify((ky.*M_y1)*(M_x2.*kx));
% c12=simplify(a12+b12);
% d12=simplify((kx.*M_x1)*(M_x2.*kx));
% e12=simplify((ky.*M_y1)*(M_y2.*ky));
% f12=simplify((kz.*M_z1)*(M_z2.*kz));
% p12=simplify((kx.*M_x1)*(M_z2.*kz));
% r12=simplify((kz.*M_z1)*(M_x2.*kx));
% s12=simplify((kz.*M_z1)*(M_y2.*ky));
% q12=simplify((ky.*M_y1)*(M_z2.*kz));
% 
% M12=c12+d12+e12+f12+p12+r12+s12+q12;
% disp('Un-Ordered Lowdin interaction matrices')
% disp(M12)

% a31=simplify((kx.*M_x3)*(M_y1.*ky));
% b31=simplify((ky.*M_y3)*(M_x1.*kx));
% c31=simplify(a31+b31);
% d31=simplify((kx.*M_x3)*(M_x1.*kx));
% e31=simplify((ky.*M_y3)*(M_y1.*ky));
% f31=simplify((kz.*M_z3)*(M_z1.*kz));
% p31=simplify((kx.*M_x3)*(M_z1.*kz));
% r31=simplify((kz.*M_z3)*(M_x1.*kx));
% s31=simplify((kz.*M_z3)*(M_y1.*ky));
% q31=simplify((ky.*M_y3)*(M_z1.*kz));
% 
% M31=c31+d31+e31+f31+p31+r31+s31+q31;
% disp('Un-Ordered Lowdin interaction matrices')
% disp(M31)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_x=sym([(sqrt(3/24)*sqrt(1/2)), (0), (sqrt(1/168)*-sqrt(1/2)), 0, 0, 0, 0, 0;
%          (0), (sqrt(15/168)*sqrt(1/2)), (0), (sqrt(3/168)*-sqrt(1/2)), 0, 0, 0, 0;
%          (0), (0), (sqrt(5/84)*sqrt(1/2)), (0), (sqrt(3/84)*-sqrt(1/2)), 0, 0, 0;
%          (0), (0), (0), (sqrt(3/84)*sqrt(1/2)), (0), (sqrt(5/84)*-sqrt(1/2)), 0, 0;
%          (0), (0), (0), (0), (sqrt(3/168)*sqrt(1/2)), (0), (sqrt(15/168)*-sqrt(1/2)), 0;
%          (0), (0), (0), (0), (0),(sqrt(1/168)*sqrt(1/2)), (0), (sqrt(3/24)*-sqrt(1/2));]);
%    
% N_y=sym([(1i*sqrt(3/24)*sqrt(1/2)), (0), (1i*sqrt(1/168)*sqrt(1/2)), 0, 0, 0, 0, 0;
%        (0), (1i*sqrt(15/168)*sqrt(1/2)), (0), (1i*sqrt(3/168)*sqrt(1/2)), 0, 0, 0, 0;
%        (0), (0), (1i*sqrt(5/84)*sqrt(1/2)), (0), (1i*sqrt(3/84)*sqrt(1/2)), 0, 0, 0;
%        (0), (0), (0), (1i*sqrt(3/84)*sqrt(1/2)), (0), (1i*sqrt(5/84)*sqrt(1/2)), 0, 0;
%        (0), (0), (0), (0), (1i*sqrt(3/168)*sqrt(1/2)), (0), (1i*sqrt(15/168)*sqrt(1/2)), 0;
%        (0), (0), (0), (0), (0),(1i*sqrt(1/168)*sqrt(1/2)), (0), (1i*sqrt(3/24)*sqrt(1/2));]);
% 
% N_z=sym([(0), (-sqrt(3/84)), (0), 0, 0, 0, 0, 0;
%        (0), (0), (-sqrt(5/84)), (0), 0, 0, 0, 0;
%        (0), (0), (0), (-sqrt(3/42)), (0), 0, 0, 0;
%        (0), (0), (0), (0), (-sqrt(3/42)), (0), 0, 0;
%        (0), (0), (0), (0), (0), (-sqrt(5/84)), (0), 0;
%        (0), (0), (0), (0), (0),(0), (-sqrt(3/84)), (0);]);
% Nx=-N_x';
% Ny=-N_y';
% Nz=-N_z'; 
%  
% disp('<Gamma_8p_5_2|Gamma_6m_7_2>')
% M_x1=simplify((R_g8m_fj5_2*(-N_x))*R_g6m_fj7_2');
% M_y1=simplify((R_g8m_fj5_2*(-N_y))*R_g6m_fj7_2');
% M_z1=simplify((R_g8m_fj5_2*(-N_z))*R_g6m_fj7_2'); 
% disp(M_x1)
% disp(M_y1)
% disp(M_z1) 
% disp('<Gamma_8p_5_2|Gamma_7p_5_2>')
% M_x2=simplify((R_g6m_fj7_2*Nx)*R_g7m_fj5_2);
% M_x3=simplify((R_g6m_fj7_2*Nx));
% M_y2=simplify((R_g6m_fj7_2*Ny)*R_g7m_fj5_2);
% M_z2=simplify((R_g6m_fj7_2*Nz)*R_g7m_fj5_2); 
% disp(M_x2)
% disp(M_x3)
% disp(M_y2)
% disp(M_z2) 
% 
% a12=simplify((kx.*M_x1)*(M_y2.*ky));
% b12=simplify((ky.*M_y1)*(M_x2.*kx));
% c12=simplify(a12+b12);
% d12=simplify((kx.*M_x1)*(M_x2.*kx));
% e12=simplify((ky.*M_y1)*(M_y2.*ky));
% f12=simplify((kz.*M_z1)*(M_z2.*kz));
% p12=simplify((kx.*M_x1)*(M_z2.*kz));
% r12=simplify((kz.*M_z1)*(M_x2.*kx));
% s12=simplify((kz.*M_z1)*(M_y2.*ky));
% q12=simplify((ky.*M_y1)*(M_z2.*kz));
% 
% M12=c12+d12+e12+f12+p12+r12+s12+q12;
% disp('Un-Ordered Lowdin interaction matrices')
% disp(M12)
