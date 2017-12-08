%%%%Eigenstates and Eigenfunctions of electrons and holes in a 
%%%%Single/Multiple Quantum Well%%%%%%%

%%%% Author: Warren Elder

%%%% This is a simple *script* file that allows the user to imput the
%%%% parameters for electrons and holes in a QW heterostructure and plots 
%%%% the solutions. The preceding files solve the exact Envelope function:
%%%%(-hbar^2/2m*)(d^2/dz^2 + k||^2)F(z)u(x,y) + V(z)u(z)u(x,y)=E.F(z)u(x,y)
%%%% Numerically using a finite difference method on a variable mesh.

clear all
clc
setup_time=0;
eigs_time=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sample_no = 'SiGe_QW_double';
%%%% Material system codes are as follows;
%%%% Si_{1-x}Ge_x,              the code is SiGe**.

Bands = 6;
substrate = 'SiGe**';
orientation = '100';
Comp_vsub = 0.0;
n_layer=5;
Material_layer1 = 'SiGe**';
Material_layer2 = 'SiGe**';
Material_layer3 = 'SiGe**';
Material_layer4 = 'SiGe**';
Material_layer5 = 'SiGe**';
Material =[ Material_layer1 
            Material_layer2 
            Material_layer3
            Material_layer4
            Material_layer5];

layer_comp=[0.91 0.91 1 0.91 0.91];
layer_npts=[100 72 100 72 100];
layer_thck=[100 72 100 72 100];
layer_cntr=[0 0 1 0 0]; % Layer/Layers of interest. If there are multiple 
                        % layers of interest, interest pointer must be 
                        % seperated by a layer.
layer_num=1:n_layer;
E_field = 0.0;          % Applied electric field [V/um]
Temp=300;               % Temperature [K]
n_Esub = 10;          	% Number of subbands
k_range= 0.06;          % Range of k_space covered
n_kpts = 3;          	% Number of k points excluding zone center
flag_1=1;               % Choose whether to save eigenvectors (0) or vector products (1)
infboundary =10;        % Choose boundary condition at layer end
%%%% Eigs parameters
options.disp=0;         % Choose 0 or 1 to not display/display eigs Ritz values
%options.p=60;         	% Choose number of Lanczos vectors
%options.tol=10E-36;       % Choose convergence tolerance. Default is{eps}~= 10E-18
%load('T_vector.mat');
%options.v0=V0;         % Choose starting vector N-by-1 vector | randomly generated
         
%%%% This section will be inputted via .txt file at a future point %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.1 Set up scaling constant
% Define fundamental constants in S.I units. All energies are in eV from
% CODATA 2006 at http://physics.nist.gov/cuu/Constants/index.html
hbar=1.054571628E-34;               %Planck's constant [Js]
m0=9.10938215E-31;                  %Free electron mass[Kg]
e=1.602176487E-19;                  %Electron charge [C]
L=1E-10;                            %Length Scale [Angstrom]
Const=(hbar*hbar)/(2*m0*e*L*L);     %Scaling constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.2 Set up substrate/virtual substrate lattice parameter

switch substrate
    case 'SiGe**'
        a0 = 5.431 + 0.1992*Comp_vsub +0.02733*Comp_vsub*Comp_vsub;
        Eg_sub = 4.2 - 3.31*Comp_vsub;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.3 Set up library and pre-allocate memory

L=repmat(struct('layer_num', 0,          ... %number of the layer
                'thick', 0,              ... %thickness of the layer
                'tot_thick', 0,          ... %total thickness of the region
                'layer_npts', 0,         ... %number of mesh points
                'tot_mesh_pts', 0,       ... %total mesh size
                'layer_mesh', 0,         ... %layer mesh sise
                'layer_cntr', 0,         ... %Layer of interest
                'comp', 0,               ... %layer composition
                'L1_vv', 0,              ... %first luttinger parameter
                'L2_vv', 0,              ... %second luttinger parameter
                'L3_vv', 0,              ... %third luttinger parameter
                'L2_vs', 0,              ... %second luttinger parameter
                'L3_vs', 0,              ... %third luttinger parameter
                'L1_ss', 0,              ... %first luttinger parameter
                'Me_g7m_cc', 0,          ... %Effective mass of electron
                'Me_g6m_cc', 0,          ... %Effective mass of electron
                'Z_g6m_vv', 0,          ... %Zeta paramter to G6- CB state
                'Z_g7m_vv', 0,          ... %Zeta paramter to G7- CB state
                'Z_g8m1_vv', 0,          ... %Zeta paramter to G8-1 CB state
                'Z_g8m2_vv', 0,          ... %Zeta paramter to G8-2 CB state
                'Z_g8m3_vv', 0,          ... %Zeta paramter to G8-3 CB state
                'Z_g7m_vs', 0,          ... %Zeta paramter to G7- CB state
                'Z_g8m1_vs', 0,          ... %Zeta paramter to G8-1 CB state
                'Z_g8m2_vs', 0,          ... %Zeta paramter to G8-2 CB state
                'Z_g7m_ss', 0,          ... %Zeta paramter to G7- CB state
                'Z_g8m1_ss', 0,          ... %Zeta paramter to G8-1 CB state
                'Z_g6p_ccg6m', 0,          ... %Zeta paramter to G6- VB state
                'Z_g8p_ccg6m', 0,          ... %Zeta paramter to G8- VB state             
                'Z_g7p_ccg7m', 0,          ... %Zeta paramter to G7- VB state
                'Z_g8p_ccg7m', 0,          ... %Zeta paramter to G8- VB state              
                'Eg_unstrained', 0,      ... %Energy gap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                'Delta', 0,              ... %spin orbit splitting
                'a',0,                   ... %Lattice parameter     
                'C11', 0,                ... %Elastic constant 'C11'
                'C12', 0,                ... %Elastic constant 'C12'                   
                'a_gap', 0,              ... %Deformation potential ('ac'-'av')
                'b', 0,                  ... %'b' Deformation potential
                'offset_vb', 0,          ... %Valence band offset
                'offset_cb', 0,          ... %Conduction band offset
                'Luttbarp', 0,            ... %Luttbar
                'Mup', 0,                 ... %Mu
                'deltap', 0,              ... %delta
                'Pip', 0,                 ... %Pi
                'Sigmap', 0,              ... %Sigma
                'exx', 0,                ... %exx
                'eyy', 0,                ... %eyy
                'ezz', 0,                ... %ezz
                'Pe_hydro', 0,           ... %Pe_hydro
                'Qe_shear', 0            ... %Qe_shear
                                            ), 1, n_layer);

% L=repmat(struct('layer_num', 0,          ... %number of the layer
%                 'thick', 0,              ... %thickness of the layer
%                 'tot_thick', 0,          ... %total thickness of the region
%                 'layer_npts', 0,         ... %number of mesh points
%                 'tot_mesh_pts', 0,       ... %total mesh size
%                 'layer_mesh', 0,         ... %layer mesh sise
%                 'layer_cntr', 0,         ... %Layer of interest
%                 'comp', 0,               ... %layer composition
%                 'Lutt1p', 0,              ... %first luttinger parameter
%                 'Lutt2p', 0,              ... %second luttinger parameter
%                 'Lutt3p', 0,              ... %third luttinger parameter
%                 'Mass_e', 0,             ... %Effective mass of electron
%                 'Eg_unstrained', 0,      ... %Energy gap
%                 'Delta', 0,              ... %spin orbit splitting
%                 'a',0,                   ... %Lattice parameter     
%                 'C11', 0,                ... %Elastic constant 'C11'
%                 'C12', 0,                ... %Elastic constant 'C12'                   
%                 'a_gap', 0,              ... %Deformation potential ('ac'-'av')
%                 'b', 0,                  ... %'b' Deformation potential
%                 'offset_vb', 0,          ... %Valence band offset
%                 'offset_cb', 0,          ... %Conduction band offset
%                 'Luttbarp', 0,            ... %Luttbar
%                 'Mup', 0,                 ... %Mu
%                 'deltap', 0,              ... %delta
%                 'Pip', 0,                 ... %Pi
%                 'Sigmap', 0,              ... %Sigma
%                 'exx', 0,                ... %exx
%                 'eyy', 0,                ... %eyy
%                 'ezz', 0,                ... %ezz
%                 'Pe_hydro', 0,           ... %Pe_hydro
%                 'Qe_shear', 0            ... %Qe_shear
%                                             ), 1, n_layer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.4 Call material parameters, set material mesh, set data file
%%%%  and pre-allocate memory
for n_1 = 1 : n_layer
    L(n_1).layer_num=layer_num(n_1); 
    L(n_1).layer_npts=layer_npts(n_1);          
    L(n_1).comp=layer_comp(n_1);
    tot_mesh_pts = sum(layer_npts); 
    tot_thick =sum(layer_thck); 
    L(n_1).layer_mesh = layer_thck(n_1)./layer_npts(n_1);
    L(n_1).layer_cntr = layer_cntr(n_1);
    L(n_1).layer_thck = layer_thck(n_1);
    L(n_1).h=L(n_1).layer_mesh;
   
    switch Material(n_1,:)

        case 'SiGe**'
        L(n_1).L1vv = (5.3484 - 0.58177*layer_comp(n_1) + 16.84*(layer_comp(n_1))^2 - 30.598*(layer_comp(n_1))^3 + 23.156*(layer_comp(n_1))^4) ; 
        L(n_1).L2vv = (0.33578 + 0.18899*layer_comp(n_1) - 8.4043*layer_comp(n_1)^2 +15.3*layer_comp(n_1)^3 - 11.578*layer_comp(n_1)^4);
        L(n_1).L3vv = (1.604 - 0.66215*layer_comp(n_1) + 8.3742*layer_comp(n_1)^2 - 15.302*layer_comp(n_1)^3 + 11.578*layer_comp(n_1)^4);
        L(n_1).L2vs = (-0.46294 + 0.12183*layer_comp(n_1) - 6.0735*layer_comp(n_1)^2 + 10.706*layer_comp(n_1)^3 - 8.3918*layer_comp(n_1)^4);
        L(n_1).L3vs = (-0.57732 + 0.072207*layer_comp(n_1) - 6.0787*layer_comp(n_1)^2 + 10.706*layer_comp(n_1)^3 - 8.3818*layer_comp(n_1)^4);
        L(n_1).L1ss = (4.007 + 1.4048*layer_comp(n_1) + 7.5549*layer_comp(n_1)^2 - 12.234*layer_comp(n_1)^3 + 10.411*layer_comp(n_1)^4) ; 
        L(n_1).Me_g7m_cc = (1/(0.528-0.49*layer_comp(n_1))); %%
        L(n_1).Me_g6m_cc = (1/(0.528-0.49*layer_comp(n_1))); %%
        L(n_1).Z_g6m_vv  = (1.3399 - 0.11*layer_comp(n_1)); 
        L(n_1).Z_g7m_vv  = (0.38018 - 0.51295*layer_comp(n_1)  - 8.3785*layer_comp(n_1)^2  + 15.302*layer_comp(n_1)^3  - 11.578*layer_comp(n_1)^4); 
        L(n_1).Z_g8m1_vv = (0.66011 + 0.23647*layer_comp(n_1) + 0.033312*layer_comp(n_1)^2 );
        L(n_1).Z_g8m2_vv = (0.058108 + 0.019301*layer_comp(n_1) + 0.002582*layer_comp(n_1)^2 );
        L(n_1).Z_g8m3_vv = (0.19585 + 0.067606*layer_comp(n_1) + 0.0092742*layer_comp(n_1)^2 );
        L(n_1).Z_g7m_vs  = (-0.51606 + 0.20568*layer_comp(n_1) - 8.5877*layer_comp(n_1)^2 + 15.14*layer_comp(n_1)^3 - 11.868*layer_comp(n_1)^4);
        L(n_1).Z_g8m1_vs = (-0.73876 - 0.24079*layer_comp(n_1) - 0.018525*layer_comp(n_1)^2 );
        L(n_1).Z_g8m2_vs = (0.21919 + 0.068556*layer_comp(n_1) + 0.0049844*layer_comp(n_1)^2 );
        L(n_1).Z_g7m_ss  = (0.69928 + 0.44415*layer_comp(n_1) - 12.219*layer_comp(n_1)^2 + 10.411*layer_comp(n_1)^3 );
        L(n_1).Z_g8m1_ss = (0.82693 + 0.24016*layer_comp(n_1) + 0.0055501*layer_comp(n_1)^2 - 0.0037868*layer_comp(n_1)^3 + 0.000034001*layer_comp(n_1)^4 );
        L(n_1).Z_g6p_ccg6m = (0);
        L(n_1).Z_g8p_ccg6m = (0);          
        L(n_1).Z_g7p_ccg7m = (0);
        L(n_1).Z_g8p_ccg7m = (0);   
        L(n_1).E_G8p = (0);
        L(n_1).E_G7p = (-0.044 - 0.2*layer_comp(n_1) - 0.052*layer_comp(n_1)^2);
        L(n_1).E_G7m = (4.15 - 3.26*layer_comp(n_1));
        L(n_1).E_G6m = (3.302 - 0.379*layer_comp(n_1));
        L(n_1).E_G8m = (3.335 - 0.222*layer_comp(n_1));
        %L(n_1).Lutt2p = (0.39 + 3.85*layer_comp(n_1)); 
        %L(n_1).Lutt3p = (1.44 + 4.13*layer_comp(n_1)); 
        %L(n_1).Mass_e = (1/(0.528-0.49*layer_comp(n_1)));
        if layer_comp(n_1)<=0.12 
        L(n_1).Eg_unstrained =3.4 -0.24*layer_comp(n_1);
        elseif layer_comp(n_1)>0.12
        L(n_1).Eg_unstrained =4.2 -3.31*layer_comp(n_1);    
        end
        %L(n_1).Delta = 0.044 + 0.245*layer_comp(n_1); 
        L(n_1).a = 5.431 + 0.1992*layer_comp(n_1)+0.02733*layer_comp(n_1)*layer_comp(n_1);
        L(n_1).C11 = 16.577 - 3.724*layer_comp(n_1);
        L(n_1).C12 = 6.393 - 1.565*layer_comp(n_1);
        L(n_1).C44 = 7.962 - 6.680*layer_comp(n_1);
        L(n_1).ac  = -10.39 - 0.02*layer_comp(n_1);
        L(n_1).av  = 1.8 - 0.56*layer_comp(n_1);
        L(n_1).a_gap = -12.19 + 0.54*layer_comp(n_1);
        L(n_1).b = -2.10 - 0.76*layer_comp(n_1);
        L(n_1).d = -4.85 - 0.43*layer_comp(n_1);
        L(n_1).offset_vb = 0.32;
        L(n_1).offset_cb = 0.68;
        %L(n_1).n_r = 3.45 + 0.55*layer_comp(n_1);
        L(n_1).E_p = 26.92 - 1.43*layer_comp(n_1);
        %L(n_1).Varshni_alpha=(0.5367+0.1475*layer_comp(n_1))/1000;
        %L(n_1).Varshni_beta=(725.8 -327.8*layer_comp(n_1))/1000;
    end 
end

n_vector=tot_mesh_pts*6;
order_ptr=zeros(n_Esub,1);
Tn_kpts=n_kpts*(n_kpts+1)/2;
Cn_kpts=0;
I1=diag(ones(6,1));
I2=diag(zeros(6,1));
infboundary_e1=infboundary;
infboundary_e2=0;
infboundary_h1=-infboundary*I1;
infboundary_h2=-infboundary*I2;
E_valelec1=ones(n_Esub,n_kpts);
E_valholes1=ones(n_Esub,n_kpts);
E_valCB{n_kpts} =zeros(n_Esub,n_kpts);
E_valVB{n_kpts} =zeros(n_Esub,n_kpts);
E_H = zeros(n_kpts,n_kpts);   
E_E = zeros(n_kpts,n_kpts);
V_hhu=zeros((tot_mesh_pts+1),n_Esub);
V_hhd=zeros((tot_mesh_pts+1),n_Esub);
V_lhu=zeros((tot_mesh_pts+1),n_Esub);
V_lhd=zeros((tot_mesh_pts+1),n_Esub);
V_sou=zeros((tot_mesh_pts+1),n_Esub);
V_sod=zeros((tot_mesh_pts+1),n_Esub);
%%%% Defines the calculated Data Structure
if flag_1==0
Data=repmat(struct( 'kx', 0,                   ... %x comp of wave vector
                    'ky', 0,                   ... %y comp of wave vector
                    'E_e', zeros(1,n_Esub),    ... %eigen energy electrons
                    'E_h', zeros(1,n_Esub),    ... %eigen energy holes
                    'V_hhu', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, HH comp up
                    'V_hhd', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, HH comp down
                    'V_lhu', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, LH comp up
                    'V_lhd', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, LH comp down
                    'V_sou', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, SO comp up
                    'V_sod', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, SO comp down
                    'V_ee', zeros(tot_mesh_pts+1, n_Esub)     ... %CB envelop fn
                            ), 1, Tn_kpts);
elseif flag_1==1
Data=repmat(struct( 'kx', 0,                        ... %x comp of wave vector
                    'ky', 0,                        ... %y comp of wave vector
                    'E_e', zeros(1,n_Esub),         ... %eigen energy electrons
                    'E_h', zeros(1,n_Esub),         ... %eigen energy holes
                    'E_cv', zeros(n_Esub,n_Esub),   ... % transition energy  
                    'chhu', zeros(n_Esub,1),        ... % coefficient of hhu 
                    'clhu', zeros(n_Esub,1),        ... % coefficient of lhu 
                    'csou', zeros(n_Esub,1),        ... % coefficient of sou 
                    'chhd', zeros(n_Esub,1),        ... % coefficient of hhd 
                    'clhd', zeros(n_Esub,1),        ... % coefficient of lhd
                    'csod', zeros(n_Esub,1),        ... % coefficient of sod 
                    'Fc_Fv', zeros(n_Esub)          ... % matrix element <electron|hole> envelope fn.
                            ), 1, Tn_kpts);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.1 Calculate the Eigenenergies and Eigenfunctions at different 
%%%% K points from X to L crystallagraphic directions
%%%% Implement waitbar
w = waitbar(0,'Please wait for me to fill up...');

for y = 0:1:n_kpts-1
for x = 0:1:n_kpts-1     
if (x>=y)
tic;
%%%% Set the in-plane dispersion and scale with the lattice constant of
%%%% substrate
Cn_kpts=Cn_kpts+1;
kx = x*(k_range/(n_kpts-1));
ky = y*(k_range/(n_kpts-1));
kpar = (2*pi/a0)*sqrt(kx*kx + ky*ky);
k_p = (2*pi/a0)*(kx + 1i*ky);
k_m = (2*pi/a0)*(kx - 1i*ky);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.2  Set symmetry coupling parameters arrising through Gamma 15 
%%%% (p-like) valence band coupling to Gamma 1 (s-like), Gamma 15 (p-like),
%%%% Gamma 12 (d-like) and Gamma 25 (f-like) conduction bands. For double 
%%%% group parameters are subsequentely adapted to treat Gamma_6, Gamma7 
%%%% Gamma8 interactions. Also set Pikus and Bir Strain Parameters
for n_2 =1:n_layer
    
% L(n_2).Luttbarp = (0.5*(L(n_2).Lutt3p + L(n_2).Lutt2p));    
% L(n_2).Mup      = (0.5*(L(n_2).Lutt3p - L(n_2).Lutt2p));
% L(n_2).deltap   = ((1/9)*(1 + L(n_2).Lutt1p + L(n_2).Lutt2p - 3*(L(n_2).Lutt3p)));
% L(n_2).Pip      = (L(n_2).Mup + (3/2)*L(n_2).deltap);
% L(n_2).Sigmap   = (L(n_2).Luttbarp - 0.5*(L(n_2).deltap));
% L(n_2).Theta    = L(n_2).Sigmap - L(n_2).deltap;
L(n_2).epar     = (a0-L(n_2).a)/L(n_2).a;
L(n_2).eperp    = -2*((L(n_2).C12)/(L(n_2).C11))*(a0-L(n_2).a)/L(n_2).a;
L(n_2).exx      = L(n_2).epar;
L(n_2).eyy      = L(n_2).epar;
L(n_2).ezz      = L(n_2).eperp;
L(n_2).exy      = 0;
L(n_2).eyz      = 0;
L(n_2).ezx      = 0;
L(n_2).Pe_hydro = L(n_2).a_gap*(L(n_2).exx + L(n_2).eyy + L(n_2).ezz);
L(n_2).Qe_shear = ((L(n_2).b)/2)*(2*L(n_2).ezz - L(n_2).exx - L(n_2).eyy);
L(n_2).Re_shear = sqrt(3)*((L(n_2).b)/2)*(L(n_2).exx - L(n_2).eyy) - 1i*L(n_2).d*L(n_2).exy;
L(n_2).Se_shear = - L(n_2).d*(L(n_2).ezx - 1i*L(n_2).eyz);

%%%% Set all mass parameters for holes to negative values
L(n_2).Lutt1    = -1*(L(n_2).Lutt1p); 
L(n_2).Lutt2    = -1*(L(n_2).Lutt2p);
L(n_2).Lutt3    = -1*(L(n_2).Lutt3p);
L(n_2).Luttbar  = -1*(L(n_2).Luttbarp);    
L(n_2).Mu       = -1*(L(n_2).Mup);
L(n_2).delta    = -1*(L(n_2).deltap);
L(n_2).Pi       = -1*(L(n_2).Pip);
L(n_2).Sigma    = -1*(L(n_2).Sigmap);
L(n_2).Theta    = -1*(L(n_2).Theta);

%%%% Defines Optical coupling between Gamma 15 and Gamma 1 bands and the
%%%% refractive index
u_2 = layer_cntr(1,n_2);
if u_2 > 0
E_p = L(n_2).E_p;
n_r = L(n_2).n_r;
end

end    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Should check this code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.3 Set Conduction and Valence band profiles

%%%%This takes the band profile and adds in the hydrostatic
%%%%and shear components.
for n_3 = 1 : n_layer    
L(n_3).Eg  =  L(n_3).Eg_unstrained + L(n_3).Pe_hydro + L(n_3).Qe_shear;
L(n_3).E_hole = L(n_3).offset_vb*(-L(n_3).Eg - (-Eg_sub));
L(n_3).E_elec = L(n_3).offset_cb*(L(n_3).Eg  - (Eg_sub) - L(n_3).Qe_shear);
if L(n_3).layer_cntr == 1    
Zeroing_1 = 0 - (L(n_3).E_hole);  
Zeroing_2 = (0 - (L(n_3).E_elec)) + L(n_3).Eg;
end
end     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Should check this code

%%%% This calculates the effect of the electric field on the band edge and
%%%% zeros all band edge to the HH
M_pt=0; T_pt=0;
for u_1 = 1 : n_layer
  L(u_1).E_hh = L(u_1).E_hole + Zeroing_1 + (L(u_1).Qe_shear);
  L(u_1).E_lh = L(u_1).E_hole + Zeroing_1 + (L(u_1).Qe_shear);
  L(u_1).E_so = L(u_1).E_hole + Zeroing_1 + (L(u_1).Qe_shear);
  L(u_1).E_el = L(u_1).E_elec + Zeroing_2;
  
%%%% This calculates the effect of the electric field on the band edge
Z_center=(tot_mesh_pts)/2;
mesh_thck=(L(u_1).layer_thck/L(u_1).layer_npts);

for u_2 = 1:(L(u_1).layer_npts)
Counter1= M_pt+ u_2;
Counter2(Counter1)= T_pt+(u_2-1)*mesh_thck;
V_E(Counter1) =E_field*(Counter2(Counter1)-Z_center)*(1e-4);
end

M_pt=M_pt+L(u_1).layer_npts;
T_pt=T_pt+L(u_1).layer_thck;
end
V_EfieldE=[V_E,0];
V_EfieldH=[V_E,0];
V_EfieldH=[V_EfieldH; V_EfieldH; V_EfieldH; V_EfieldH; V_EfieldH; V_EfieldH];
V_EfieldH=(diag(reshape(V_EfieldH,1,(6*(tot_mesh_pts+1)))));

if Cn_kpts == 1
B_P.V_Efield=V_E;
for u_3 = 1:n_layer
u_4 = layer_cntr(1,u_3);
if u_4>0
Layermeshpts=L(u_3).layer_npts;
B_P(u_3).Band_profile_1 =  [(L(u_3).E_hh - L(u_3).Qe_shear)   (0)                              (0)
                            (0)                               (L(u_3).E_lh + L(u_3).Qe_shear)  (-sqrt(2)*L(u_3).Qe_shear)
                            (0)                               (-sqrt(2)*L(u_3).Qe_shear)       (L(u_3).E_so - L(u_3).Delta)];                           
B_P(u_3).Band_profile_2 =  (eig(B_P(u_3).Band_profile_1));                                                
B_P(u_3).E_heavyhole= B_P(u_3).Band_profile_2(3)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole= B_P(u_3).Band_profile_2(2)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit= B_P(u_3).Band_profile_2(1)*ones(1,Layermeshpts);
else
Layermeshpts=L(u_3).layer_npts;
B_P(u_3).Band_profile_1 =  [(L(u_3).E_hh - L(u_3).Qe_shear)   (0)                              (0)
                            (0)                               (L(u_3).E_lh + L(u_3).Qe_shear)  (-sqrt(2)*L(u_3).Qe_shear)
                            (0)                               (-sqrt(2)*L(u_3).Qe_shear)       (L(u_3).E_so - L(u_3).Delta)];                            
B_P(u_3).Band_profile_2 =  (eig(B_P(u_3).Band_profile_1));                                              
B_P(u_3).E_heavyhole= B_P(u_3).Band_profile_2(2)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole= B_P(u_3).Band_profile_2(3)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit= B_P(u_3).Band_profile_2(1)*ones(1,Layermeshpts);
end
end
end    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Should check this code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.4 This loop calculate the Block Tridiagonal Matrix Hamiltonian 
%%%% elements. We are using the forward-backward finite difference 
%%%% approximation to obtain discretitized values of the Hamiltonioan at 
%%%% each of the N mesh points. These values are then re-cast into a cell 
%%%% array forming a block tridiagonal matrix consisting of cells for 
%%%% diagonal, sub and super diagonal at each mesh point. Set the 
%%%% Hamiltonian elements P, Q, R, S, Epsilon and C where h1=-id/dz and 
%%%% h2=d^2/dz^2. Below are the paramaters before discretisation by finite
%%%% element theorum. These are taken from Foreman 93.
 
%%%%P = Pe + Const*((Lutt1*kpar.*kpar + h2*Lutt1));
%%%%Q = Qe + Const*((Lutt2*kpar.*kpar - h2*2*Lutt2));  
%%%%S = Const*(2*sqrt(3)*kpar.*((Sigma - delta)*h1 + h1*(Pi)));
%%%%Epsilon = Const*(2*sqrt(3)*kpar.*(((1/3)*(Sigma - delta) +
%%%%         (2/3)*(Pi))*h1 + h1*((2/3)*(Sigma - delta) + (1/3)*(Pi))));
%%%%R = Const*(-(sqrt(3))*Luttphi*kpar.*kpar);
%%%%C = Const*(2*kpar.*((h1*(Sigma - delta - Pi) - (Sigma - delta -
%%%%    Pi)*h1)));

for n_4 = 1 : n_layer
         
    if (n_4 == 1) % Set up first interface matrix at start of structure
    %%%% Stage 1 Center Fj=i points 
    L(n_4).E_elI= infboundary;
    L(n_4).E_hhI=-infboundary;
    L(n_4).E_lhI=-infboundary;
    L(n_4).E_soI=-infboundary;        
    L(n_4).DeltaI = 0;
    L(n_4).Qe_shearI = 0;
    L(n_4).M_ecI= 0;
    L(n_4).PcI = 0;      
    L(n_4).QcI = 0;      
    L(n_4).RcI = 0;
    L(n_4).RcI_dag = 0;
    %%%% Stage 2 Forward F{j=i+1) points
    L(n_4).M_epI= 0;
    L(n_4).PpI = 0;
    L(n_4).QpI = 0;     
    L(n_4).SpI = 0;
    L(n_4).SpI_dag = 0;  
    L(n_4).EpsilonpI = 0;
    L(n_4).EpsilonpI_dag = 0;
    L(n_4).CpI = 0;
    L(n_4).CpI_dag = 0;
    %%%% Stage 3 Backward F{j=i-1} points
    L(n_4).M_emI= 0;        
    L(n_4).PmI = 0;
    L(n_4).QmI = 0;   
    L(n_4).SmI = 0;
    L(n_4).SmI_dag = 0;
    L(n_4).EpsilonmI = 0;
    L(n_4).EpsilonmI_dag = 0;
    L(n_4).CmI = 0;
    L(n_4).CmI_dag = 0;
    %%%%% Stage 4 Center F{j=i-1} Interface LEFT
    L(n_4).E_elIL=0;
    L(n_4).E_hhIL=0;
    L(n_4).E_lhIL=0;
    L(n_4).E_soIL=0;
    L(n_4).PcIL = 0;  
    L(n_4).QcIL = 0;     
    L(n_4).DeltaIL = 0;
    L(n_4).Qe_shearIL = 0; 
    L(n_4).RcIL = 0;
    L(n_4).RcIL_dag = 0;
    %%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
    L(n_4).E_elIR=0;
    L(n_4).E_hhIR=0;
    L(n_4).E_lhIR=0;
    L(n_4).E_soIR=0;
    L(n_4).PcIR = 0;  
    L(n_4).QcIR = 0;     
    L(n_4).DeltaIR = 0;
    L(n_4).Qe_shearIR = 0; 
    L(n_4).RcIR = 0;
    L(n_4).RcIR_dag = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Set up other interface matrix
    
    if (L(n_4).layer_cntr == 1) || (L(n_4-1).layer_cntr == 1) 
    %%%% These are interface matrices of interest as they coincide with the
    %%%% abrupt change in the material parameters.
    %%%% Stage 1 Center F{j=i} points    
    L(n_4).E_elI=((L(n_4).E_el + L(n_4-1).E_el)/2);
    L(n_4).E_hhI=((L(n_4).E_hh + L(n_4-1).E_hh)/2);
    L(n_4).E_lhI=((L(n_4).E_lh + L(n_4-1).E_lh)/2);
    L(n_4).E_soI=((L(n_4).E_so + L(n_4-1).E_so)/2);
    L(n_4).M_ecI= Const*(kpar*kpar*((((L(n_4).Mass_e)+(L(n_4-1).Mass_e))/2)) - (-2)*(((L(n_4).Mass_e)*L(n_4-1).h...
               + (L(n_4-1).Mass_e)*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)))); 
    
    L(n_4).PcI = Const*(kpar*kpar*((((L(n_4).Lutt1)+(L(n_4-1).Lutt1))/2)) - (-2)*(((L(n_4).Lutt1).*L(n_4-1).h...
               + (L(n_4-1).Lutt1)*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))));  
           
    L(n_4).QcI = Const*(kpar*kpar*((((L(n_4).Lutt2)+(L(n_4-1).Lutt2))/2)) - (-2)*((((-2)*L(n_4).Lutt2)*L(n_4-1).h...
               + ((-2)*L(n_4-1).Lutt2).*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)))); 
           
    L(n_4).DeltaI = ((L(n_4).Delta +L(n_4-1).Delta)/2);
    L(n_4).Qe_shearI = ((L(n_4).Qe_shear + L(n_4-1).Qe_shear)/2); 
    L(n_4).RcI = Const*sqrt(3)*((-((L(n_4).Luttbar + L(n_4-1).Luttbar)/2)*k_m*k_m)+(((L(n_4).Mu + L(n_4-1).Mu)/2)*k_p*k_p));
    L(n_4).RcI_dag = Const*sqrt(3)*((-((L(n_4).Luttbar + L(n_4-1).Luttbar)/2)*k_p*k_p)+(((L(n_4).Mu + L(n_4-1).Mu)/2)*k_m*k_m));
   
    %%%% Stage 2 Forward F{j=i+1) points
    L(n_4).M_epI= Const*(-2*(L(n_4).Mass_e*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));    
    L(n_4).PpI = Const*(-2*((L(n_4).Lutt1+((L(n_4-1).Lutt1-L(n_4).Lutt1)/4))*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));   
    L(n_4).QpI = Const*(-2*(-2*(L(n_4).Lutt2+((L(n_4-1).Lutt2-L(n_4).Lutt2)/4)).*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))); 
    L(n_4).SpI = Const*2*sqrt(3)*...
        (((((L(n_4-1).Theta)) + ((L(n_4).Theta)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (L(n_4).Pi./(L(n_4).h + L(n_4-1).h)));
    L(n_4).SpI_dag = Const*2*sqrt(3)*...
        (((((L(n_4-1).Pi)) + ((L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (L(n_4).Theta./(L(n_4).h + L(n_4-1).h))); 
    L(n_4).EpsilonpI = Const*2*sqrt(3)*(...
        ((((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))...
        + ((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
        + (((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi))/(L(n_4).h + L(n_4-1).h)));
    L(n_4).EpsilonpI_dag = Const*2*sqrt(3)*(...
        ((((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/((L(n_4).h + L(n_4-1).h)))...
        + (((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi)) ...
        + ((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)));
    L(n_4).CpI = Const*2*...
        (((((L(n_4).Theta) - L(n_4).Pi) - ((L(n_4-1).Theta) - L(n_4-1).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CpI_dag = Const*2*...
        (((((L(n_4-1).Theta) - L(n_4-1).Pi) - ((L(n_4).Theta) - L(n_4).Pi))/(2*(L(n_4).h + L(n_4-1).h))));

    %%%% Stage 3 Backward F{j=i-1} points
    L(n_4).M_emI= Const*(-2*(L(n_4-1).Mass_e.*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));        
    L(n_4).PmI  = Const*(-2*((L(n_4-1).Lutt1-((L(n_4-1).Lutt1-L(n_4).Lutt1)/4))*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));
    L(n_4).QmI  = Const*(-2*(-2*(L(n_4-1).Lutt2-((L(n_4-1).Lutt2-L(n_4).Lutt2)/4))*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));   
    L(n_4).SmI = Const*2*sqrt(3)*...
        (-((((L(n_4-1).Theta)) + ((L(n_4).Theta)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (-(L(n_4-1).Pi/(L(n_4).h + L(n_4-1).h))));
    L(n_4).SmI_dag = Const*2*sqrt(3)*...
        (-((((L(n_4-1).Pi)) + ((L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (-(L(n_4-1).Theta/(L(n_4).h + L(n_4-1).h))));  
    L(n_4).EpsilonmI = Const*2*sqrt(3)*(...
        (-(((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))...
        + ((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
        + (-(((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi))/(L(n_4).h + L(n_4-1).h))));
    L(n_4).EpsilonmI_dag = Const*2*sqrt(3)*(...
        (-((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))/((L(n_4).h + L(n_4-1).h)))...
        + (-((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi)...
        + ((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CmI = Const*2*...
        (((((L(n_4).Theta) - L(n_4).Pi) - ((L(n_4-1).Theta) - L(n_4-1).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CmI_dag = Const*2*...
        (((((L(n_4-1).Theta) - L(n_4-1).Pi) - ((L(n_4).Theta) - L(n_4).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    
    %%%%% Stage 4 Center F{j=i-1} Interface LEFT 
    L(n_4).E_elIL=(L(n_4-1).E_el);
    L(n_4).E_hhIL=(L(n_4-1).E_hh);
    L(n_4).E_lhIL=(L(n_4-1).E_lh);
    L(n_4).E_soIL=(L(n_4-1).E_so);
    L(n_4).PcIL = Const*(kpar*kpar*(L(n_4-1).Lutt1)-...
                 (-2)*((1.5*(L(n_4-1).Lutt1)+0.5*((L(n_4).Lutt1 + L(n_4-1).Lutt1)/2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));  
    L(n_4).QcIL = Const*(kpar*kpar*(L(n_4-1).Lutt2)-...
                 (-2)*((-2)*(1.5*(L(n_4-1).Lutt2)+0.5*((L(n_4).Lutt2 + L(n_4-1).Lutt2)/2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));     
    L(n_4).DeltaIL = (L(n_4-1).Delta);
    L(n_4).Qe_shearIL = (L(n_4-1).Qe_shear); 
    L(n_4).RcIL = Const*sqrt(3)*((-(L(n_4-1).Luttbar)*k_m*k_m)+((L(n_4-1).Mu)*k_p*k_p));
    L(n_4).RcIL_dag = Const*sqrt(3)*((-(L(n_4-1).Luttbar)*k_p*k_p)+((L(n_4-1).Mu)*k_m*k_m));
    
    %%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
    L(n_4).E_elIR=(L(n_4).E_el);
    L(n_4).E_hhIR=(L(n_4).E_hh);
    L(n_4).E_lhIR=(L(n_4).E_lh);
    L(n_4).E_soIR=(L(n_4).E_so);
    L(n_4).PcIR = Const*(kpar*kpar*(L(n_4).Lutt1)-...
                 (-2)*((1.5*(L(n_4).Lutt1)+0.5*((L(n_4).Lutt1 + L(n_4-1).Lutt1)/2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));  
    
    L(n_4).QcIR = Const*(kpar*kpar*(L(n_4).Lutt2)-...
                 (-2)*((-2)*(1.5*(L(n_4).Lutt2)+0.5*((L(n_4).Lutt2 + L(n_4-1).Lutt2)/2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));  
             
    L(n_4).DeltaIR = (L(n_4).Delta);
    L(n_4).Qe_shearIR = (L(n_4).Qe_shear); 
    L(n_4).RcIR = Const*sqrt(3)*((-(L(n_4).Luttbar)*k_m*k_m)+((L(n_4).Mu)*k_p*k_p));
    L(n_4).RcIR_dag = Const*sqrt(3)*((-(L(n_4).Luttbar)*k_p*k_p)+((L(n_4).Mu)*k_m*k_m));
      
    else % All other quasi-interfacial points within a layer
    
    %%%% Stage 1 Center F{j=i} points    
    L(n_4).E_elI=L(n_4).E_el;
    L(n_4).E_hhI=L(n_4).E_hh;
    L(n_4).E_lhI=L(n_4).E_lh;
    L(n_4).E_soI=L(n_4).E_so;
    L(n_4).DeltaI = L(n_4).Delta;
    L(n_4).M_ecI= Const*(kpar*kpar*(L(n_4).Mass_e) - (-2)*(((L(n_4).Mass_e)*L(n_4-1).h...
               + (L(n_4-1).Mass_e)*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))));   
    L(n_4).PcI = Const*(kpar*kpar*(L(n_4).Lutt1) - (-2)*(((L(n_4).Lutt1).*L(n_4-1).h...
               + (L(n_4-1).Lutt1)*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))));       
    L(n_4).QcI = Const*(kpar*kpar*(L(n_4).Lutt2) - (-2)*((((-2)*L(n_4).Lutt2)*L(n_4-1).h...
               + ((-2)*L(n_4-1).Lutt2).*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))));      
    L(n_4).Qe_shearI = L(n_4).Qe_shear;
    L(n_4).RcI = Const*sqrt(3)*((-L(n_4).Luttbar.*k_m.*k_m)+(L(n_4).Mu*k_p*k_p));
    L(n_4).RcI_dag = Const*sqrt(3)*((-L(n_4).Luttbar.*k_p.*k_p)+(L(n_4).Mu*k_m*k_m));
    
    %%%% Stage 2 Forward F{j=i+1) points
    L(n_4).M_epI= Const*(-2*(L(n_4).Mass_e*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));    
    L(n_4).PpI = Const*(-2*(L(n_4).Lutt1*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));   
    L(n_4).QpI = Const*(-2*(-2*L(n_4).Lutt2.*L(n_4).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h))); 
    L(n_4).SpI = Const*2*sqrt(3)*...
        (((((L(n_4-1).Theta)) + ((L(n_4).Theta)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (L(n_4).Pi./(L(n_4).h + L(n_4-1).h)));
    L(n_4).SpI_dag = Const*2*sqrt(3)*...
        (((((L(n_4-1).Pi)) + ((L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (L(n_4).Theta./(L(n_4).h + L(n_4-1).h))); 
    L(n_4).EpsilonpI = Const*2*sqrt(3)*(...
        ((((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))...
        + ((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
        + (((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi))/(L(n_4).h + L(n_4-1).h)));
    L(n_4).EpsilonpI_dag = Const*2*sqrt(3)*(...
        ((((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/((L(n_4).h + L(n_4-1).h)))...
        + (((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi)) ...
        + ((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)));
    L(n_4).CpI = Const*2*...
        (((((L(n_4).Theta) - L(n_4).Pi) - ((L(n_4-1).Theta) - L(n_4-1).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CpI_dag = Const*2*...
        (((((L(n_4-1).Theta) - L(n_4-1).Pi) - ((L(n_4).Theta) - L(n_4).Pi))/(2*(L(n_4).h + L(n_4-1).h))));

    %%%% Stage 3 Backward F{j=i-1} points
    L(n_4).M_emI= Const*(-2*(L(n_4-1).Mass_e.*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));        
    L(n_4).PmI  = Const*(-2*(L(n_4-1).Lutt1*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));
    L(n_4).QmI  = Const*(-2*(-2*L(n_4-1).Lutt2*L(n_4-1).h)/(L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)));   
    L(n_4).SmI = Const*2*sqrt(3)*...
        (-((((L(n_4-1).Theta)) + ((L(n_4).Theta)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (-(L(n_4-1).Pi/(L(n_4).h + L(n_4-1).h))));
    L(n_4).SmI_dag = Const*2*sqrt(3)*...
        (-((((L(n_4-1).Pi)) + ((L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
         + (-(L(n_4-1).Theta/(L(n_4).h + L(n_4-1).h))));  
    L(n_4).EpsilonmI = Const*2*sqrt(3)*(...
        (-(((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))...
        + ((1/3)*((L(n_4).Theta)) + (2/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h)))...
        + (-(((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi))/(L(n_4).h + L(n_4-1).h))));
    L(n_4).EpsilonmI_dag = Const*2*sqrt(3)*(...
        (-((1/3)*((L(n_4-1).Theta)) + (2/3)*(L(n_4-1).Pi))/((L(n_4).h + L(n_4-1).h)))...
        + (-((2/3)*((L(n_4-1).Theta)) + (1/3)*(L(n_4-1).Pi)...
        + ((2/3)*((L(n_4).Theta)) + (1/3)*(L(n_4).Pi)))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CmI = Const*2*...
        (((((L(n_4).Theta) - L(n_4).Pi) - ((L(n_4-1).Theta) - L(n_4-1).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    L(n_4).CmI_dag = Const*2*...
        (((((L(n_4-1).Theta) - L(n_4-1).Pi) - ((L(n_4).Theta) - L(n_4).Pi))/(2*(L(n_4).h + L(n_4-1).h))));
    
    %%%%% Stage 4 Center F{j=i-1} Interface LEFT
    L(n_4).E_elIL=(L(n_4-1).E_el);
    L(n_4).E_hhIL=(L(n_4-1).E_hh);
    L(n_4).E_lhIL=(L(n_4-1).E_lh);
    L(n_4).E_soIL=(L(n_4-1).E_so);
    L(n_4).PcIL = Const*(kpar*kpar*(L(n_4-1).Lutt1)-...
                 (-2)*((1.5*(L(n_4-1).Lutt1)+0.5*(L(n_4).Lutt1 + L(n_4-1).Lutt1))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));  
    L(n_4).QcIL = Const*(kpar*kpar*(L(n_4-1).Lutt2)-...
                 (-2)*((-2)*(1.5*(L(n_4-1).Lutt2)+0.5*(L(n_4).Lutt2 + L(n_4-1).Lutt2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));     
    L(n_4).DeltaIL = (L(n_4-1).Delta);
    L(n_4).Qe_shearIL = (L(n_4-1).Qe_shear); 
    L(n_4).RcIL = Const*sqrt(3)*((-(L(n_4-1).Luttbar)*k_m*k_m)+((L(n_4-1).Mu)*k_p*k_p));
    L(n_4).RcIL_dag = Const*sqrt(3)*((-(L(n_4-1).Luttbar)*k_p*k_p)+((L(n_4-1).Mu)*k_m*k_m));
    
    %%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
    L(n_4).E_elIR=(L(n_4).E_el);
    L(n_4).E_hhIR=(L(n_4).E_hh);
    L(n_4).E_lhIR=(L(n_4).E_lh);
    L(n_4).E_soIR=(L(n_4).E_so);
    L(n_4).PcIR = Const*(kpar*kpar*(L(n_4).Lutt1)-...
                 (-2)*((1.5*(L(n_4).Lutt1)+0.5*(L(n_4).Lutt1 + L(n_4-1).Lutt1))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));  
    L(n_4).QcIR = Const*(kpar*kpar*(L(n_4).Lutt2)-...
                 (-2)*((-2)*(1.5*(L(n_4).Lutt2)+0.5*(L(n_4).Lutt2 + L(n_4-1).Lutt2))...
                 /(L(n_4).h*(L(n_4).h+L(n_4).h))));     
    L(n_4).DeltaIR = (L(n_4).Delta);
    L(n_4).Qe_shearIR = (L(n_4).Qe_shear); 
    L(n_4).RcIR = Const*sqrt(3)*((-(L(n_4).Luttbar)*k_m*k_m)+((L(n_4).Mu)*k_p*k_p));
    L(n_4).RcIR_dag = Const*sqrt(3)*((-(L(n_4).Luttbar)*k_p*k_p)+((L(n_4).Mu)*k_m*k_m));
    
    end
    end
    
    %%%% Set up all other j points away from interface in layer
    %%%% Stage 1 Center Fj=i points
    L(n_4).E_elJ=L(n_4).E_el;
    L(n_4).E_hhJ=L(n_4).E_hh;
    L(n_4).E_lhJ=L(n_4).E_lh;
    L(n_4).E_soJ=L(n_4).E_so;
    L(n_4).DeltaJ = L(n_4).Delta;
    L(n_4).Qe_shearJ = L(n_4).Qe_shear;
    L(n_4).M_ecJ= Const*(kpar*kpar*(L(n_4).Mass_e) - ...
        (-2)*(((L(n_4).Mass_e).*L(n_4).h + (L(n_4).Mass_e)*L(n_4).h)/(L(n_4).h*L(n_4).h*(L(n_4).h+L(n_4).h))));   
    L(n_4).PcJ = Const*(kpar*kpar*(L(n_4).Lutt1) - ...
        (-2)*(((L(n_4).Lutt1).*L(n_4).h + (L(n_4).Lutt1)*L(n_4).h)/(L(n_4).h*L(n_4).h*(L(n_4).h+L(n_4).h)))); 
    L(n_4).QcJ = Const*(kpar*kpar*(L(n_4).Lutt2) - ...
        (-2)*((((-2)*L(n_4).Lutt2)*L(n_4).h + ((-2)*L(n_4).Lutt2)*L(n_4).h)/(L(n_4).h*L(n_4).h*(L(n_4).h+L(n_4).h))));
    L(n_4).RcJ = Const*sqrt(3)*((-L(n_4).Luttbar.*k_m*k_m)+(L(n_4).Mu*k_p*k_p));
    L(n_4).RcJ_dag = Const*sqrt(3)*((-L(n_4).Luttbar.*k_p*k_p)+(L(n_4).Mu*k_m*k_m));

    %%%% Stage 2 Forward F{j=i+1) points
    L(n_4).M_epJ= Const*(-2*(L(n_4).Mass_e*L(n_4).h)/(L(n_4).h*L(n_4).h*(L(n_4).h+L(n_4).h)));        
    L(n_4).PpJ = Const*(-2*(L(n_4).Lutt1.*L(n_4).h)/(L(n_4).h.*L(n_4).h.*(L(n_4).h+L(n_4).h)));
    L(n_4).QpJ = Const*(-2*(-2*L(n_4).Lutt2*L(n_4).h)/(L(n_4).h.*L(n_4).h.*(L(n_4).h+L(n_4).h)));     
    L(n_4).SpJ = Const*2*sqrt(3)*...
        (((L(n_4).Theta+L(n_4).Pi)./(L(n_4).h + L(n_4).h)));   
    L(n_4).EpsilonpJ = Const*2*sqrt(3)*(...
        ((((1/3)*(L(n_4).Theta) + (2/3)*(L(n_4).Pi)))/(L(n_4).h + L(n_4).h))...
        + (((2/3)*(L(n_4).Theta) + (1/3)*(L(n_4).Pi))/(L(n_4).h + L(n_4).h)));
    
    %%%% Stage 3 Backward F{j=i-1} points 
    L(n_4).M_emJ= Const*(-2*(L(n_4).Mass_e*L(n_4).h)/(L(n_4).h.*L(n_4).h.*(L(n_4).h+L(n_4).h)));
    L(n_4).PmJ = Const*(-2*(L(n_4).Lutt1*L(n_4).h)/(L(n_4).h.*L(n_4).h.*(L(n_4).h+L(n_4).h)));
    L(n_4).QmJ = Const*(-2*(-2*L(n_4).Lutt2*L(n_4).h)/(L(n_4).h.*L(n_4).h.*(L(n_4).h+L(n_4).h)));   
    L(n_4).SmJ = Const*2*sqrt(3)*...
        (-((L(n_4).Theta + L(n_4).Pi)./(L(n_4).h + L(n_4).h)));  
    L(n_4).EpsilonmJ = Const*2*sqrt(3)*(...
        (-(((1/3)*(L(n_4).Theta) + (2/3)*(L(n_4).Pi)))/(L(n_4).h + L(n_4).h))...
        + (-(((2/3)*(L(n_4).Theta) + (1/3)*(L(n_4).Pi))/(L(n_4).h + L(n_4).h))));
%Define the Hamiltonian that will be solved in the finite difference
%theorum for Electrons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Strain Hamiltonian

% L(n_4).Hevv_centerI =  diag([(-L(n_4).Qe_I)   
%                              (L(n_4).Qe_I)
%                              (L(n_4).Qe_I)
%                              (-L(n_4).Qe_I)]);  
%                          
% L(n_4).Hevv_centerIR =  diag([(-L(n_4).Qe_IR)   
%                              (L(n_4).Qe_IR)
%                              (L(n_4).Qe_IR)
%                              (-L(n_4).Qe_IR)]); 
%                          
% L(n_4).Hevv_centerIL =  diag([(-L(n_4).Qe_IL)   
%                              (L(n_4).Qe_IL)
%                              (L(n_4).Qe_IL)
%                              (-L(n_4).Qe_IL)]); 
%                          
% L(n_4).Hevv_centerI =  diag([(-L(n_4).Qe_J)   
%                              (L(n_4).Qe_J)
%                              (L(n_4).Qe_J)
%                              (-L(n_4).Qe_J)]);                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Band offset energies
L(n_4).Ecc_centerI  =  diag([(L(n_4).Ecc_cI)  
                             (L(n_4).Ecc_cI) ]);
                        
L(n_4).Ecc_centerIR =  diag([(L(n_4).Ecc_cIR)  
                             (L(n_4).Ecc_cIR) ]);
                        
L(n_4).Ecc_centerIL =  diag([(L(n_4).Ecc_cIL)  
                             (L(n_4).Ecc_cIL) ]);

L(n_4).Ecc_centerJ  =  diag([(L(n_4).Ecc_cJ)  
                             (L(n_4).Ecc_cJ) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
L(n_4).Evv_centerI  =  diag([(L(n_4).Evv_cI)  
                             (L(n_4).Evv_cI)
                             (L(n_4).Evv_cI)
                             (L(n_4).Evv_cI)]);
                        
L(n_4).Evv_centerIR =  diag([(L(n_4).Evv_cIR)  
                             (L(n_4).Evv_cIR)
                             (L(n_4).Evv_cIR)
                             (L(n_4).Evv_cIR)]);
                        
L(n_4).Evv_centerIL =  diag([(L(n_4).Evv_cIL)  
                             (L(n_4).Evv_cIL)
                             (L(n_4).Evv_cIL)
                             (L(n_4).Evv_cIL)]);
                         
L(n_4).Evv_centerJ  =  diag([(L(n_4).Evv_cJ)  
                             (L(n_4).Evv_cJ)
                             (L(n_4).Evv_cJ)
                             (L(n_4).Evv_cJ)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
L(n_4).Ess_centerI  =  diag([(L(n_4).Ess_cI)  
                             (L(n_4).Ess_cI) ]);

L(n_4).Ess_centerIR =  diag([(L(n_4).Ess_cIR)  
                             (L(n_4).Ess_cIR) ]);

L(n_4).Ess_centerIL =  diag([(L(n_4).Ess_cIL)  
                             (L(n_4).Ess_cIL) ]);
                        
L(n_4).Ess_centerJ  =  diag([(L(n_4).Ess_cJ)  
                             (L(n_4).Ess_cJ) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%%%% Conduction band block
%%%% Interface
L(n_4).Lcc_centerI =  [(L(n_4).Pe_cI)     (L(n_4).Ce_cI) 
                       (L(n_4).Ce_cI_dag) (L(n_4).Pe_cI_con) ];
                   
L(n_4).Lcc_centerIR =  [(L(n_4).Pe_cIR)     (L(n_4).Ce_cIR) 
                       (L(n_4).Ce_cIR_dag) (L(n_4).Pe_cIR_con) ];
                   
L(n_4).Lcc_centerIL =  [(L(n_4).Pe_cIL)     (L(n_4).Ce_cIL) 
                       (L(n_4).Ce_cIL_dag) (L(n_4).Pe_cIL_con) ];

L(n_4).Lcc_plusI   =  [(L(n_4).Pe_pI)     (L(n_4).Ce_pI) 
                       (L(n_4).Ce_pI_dag) (L(n_4).Pe_pI_con) ];

L(n_4).Lcc_minusI  =  [(L(n_4).Pe_mI)     (L(n_4).Ce_mI) 
                       (L(n_4).Ce_mI_dag) (L(n_4).Pe_mI_con) ];
%%%% LAYER                    
L(n_4).Lcc_centerJ =  [(L(n_4).Pe_cJ)     (L(n_4).Ce_cJ) 
                       (L(n_4).Ce_cJ_dag) (L(n_4).Pe_cJ_con) ];

L(n_4).Lcc_plusJ   =  [(L(n_4).Pe_pJ)     (L(n_4).Ce_pJ) 
                       (L(n_4).Ce_pJ_dag) (L(n_4).Pe_pJ_con) ];

L(n_4).Lcc_minusJ  =  [(L(n_4).Pe_mJ)     (L(n_4).Ce_mJ) 
                       (L(n_4).Ce_mJ_dag) (L(n_4).Pe_mJ_con) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
%%%% Valence band block
%%%% Interface
L(n_4).Lvv_centerI =  [(L(n_4).Pvv_cI + L(n_4).Qvv_cI) (L(n_4).Svv_cI)                 (L(n_4).Rvv_cI)                         (L(n_4).Cvv_cI_dag)
                       (L(n_4).Svv_cI_dag)             (L(n_4).Pvv_cI - L(n_4).Qvv_cI) (L(n_4).Zvv_cI)                         (L(n_4).Rvv_cI)
                       (L(n_4).Rvv_cI_dag)             (L(n_4).Zvv_cI_dag)             (L(n_4).Pvv_cI_con - L(n_4).Qvv_cI_con) (L(n_4).Svv_cI_dag_con)
                       (L(n_4).Cvv_cI)                 (L(n_4).Rvv_cI_dag)             (L(n_4).Svv_cI_con)                     (L(n_4).Pvv_cI_con - L(n_4).Qvv_cI_con) ];

L(n_4).Lvv_centerIR =  [(L(n_4).Pvv_cIR + L(n_4).Qvv_cIR) (L(n_4).Svv_cIR)                  (L(n_4).Rvv_cIR)                          (L(n_4).Cvv_cIR_dag)
                       (L(n_4).Svv_cIR_dag)               (L(n_4).Pvv_cIR - L(n_4).Qvv_cIR) (L(n_4).Zvv_cIR)                          (L(n_4).Rvv_cIR)
                       (L(n_4).Rvv_cIR_dag)               (L(n_4).Zvv_cIR_dag)              (L(n_4).Pvv_cIR_con - L(n_4).Qvv_cIR_con) (L(n_4).Svv_cIR_dag_con)
                       (L(n_4).Cvv_cIR)                   (L(n_4).Rvv_cIR_dag)              (L(n_4).Svv_cIR_con)                      (L(n_4).Pvv_cIR_con - L(n_4).Qvv_cIR_con) ];
                   
L(n_4).Lvv_centerIL =  [(L(n_4).Pvv_cIL + L(n_4).Qvv_cIL) (L(n_4).Svv_cIL)                  (L(n_4).Rvv_cIL)                          (L(n_4).Cvv_cIL_dag)
                       (L(n_4).Svv_cIL_dag)               (L(n_4).Pvv_cIL - L(n_4).Qvv_cIL) (L(n_4).Zvv_cIL)                          (L(n_4).Rvv_cIL)
                       (L(n_4).Rvv_cIL_dag)               (L(n_4).Zvv_cIL_dag)              (L(n_4).Pvv_cIL_con - L(n_4).Qvv_cIL_con) (L(n_4).Svv_cIL_dag_con)
                       (L(n_4).Cvv_cIL)                   (L(n_4).Rvv_cIL_dag)              (L(n_4).Svv_cIL_con)                      (L(n_4).Pvv_cIL_con - L(n_4).Qvv_cIL_con) ];
                                     
L(n_4).Lvv_plusI =    [(L(n_4).Pvv_pI + L(n_4).Qvv_pI) (L(n_4).Svv_pI)                 (L(n_4).Rvv_pI)                         (L(n_4).Cvv_pI_dag)
                       (L(n_4).Svv_pI_dag)             (L(n_4).Pvv_pI - L(n_4).Qvv_pI) (L(n_4).Zvv_pI)                         (L(n_4).Rvv_pI)
                       (L(n_4).Rvv_pI_dag)             (L(n_4).Zvv_pI_dag)             (L(n_4).Pvv_pI_con - L(n_4).Qvv_pI_con) (L(n_4).Svv_pI_dag_con)
                       (L(n_4).Cvv_pI)                 (L(n_4).Rvv_pI_dag)             (L(n_4).Svv_pI_con)                     (L(n_4).Pvv_pI_con - L(n_4).Qvv_pI_con) ];

L(n_4).Lvv_minusI =   [(L(n_4).Pvv_mI + L(n_4).Qvv_mI) (L(n_4).Svv_mI)                 (L(n_4).Rvv_mI)                         (L(n_4).Cvv_mI_dag)
                       (L(n_4).Svv_mI_dag)             (L(n_4).Pvv_mI - L(n_4).Qvv_mI) (L(n_4).Zvv_mI)                         (L(n_4).Rvv_mI)
                       (L(n_4).Rvv_mI_dag)             (L(n_4).Zvv_mI_dag)             (L(n_4).Pvv_mI_con - L(n_4).Qvv_mI_con) (L(n_4).Svv_mI_dag_con)
                       (L(n_4).Cvv_mI)                 (L(n_4).Rvv_mI_dag)             (L(n_4).Svv_mI_con)                     (L(n_4).Pvv_mI_con - L(n_4).Qvv_mI_con) ];
                   
%%%% LAYER                    
L(n_4).Lvv_centerJ =  [(L(n_4).Pvv_cJ + L(n_4).Qvv_cJ) (L(n_4).Svv_cJ)                 (L(n_4).Rvv_cJ)                         (L(n_4).Cvv_cJ_dag)
                       (L(n_4).Svv_cJ_dag)             (L(n_4).Pvv_cJ - L(n_4).Qvv_cJ) (L(n_4).Zvv_cJ)                         (L(n_4).Rvv_cJ)
                       (L(n_4).Rvv_cJ_dag)             (L(n_4).Zvv_cJ_dag)             (L(n_4).Pvv_cJ_con - L(n_4).Qvv_cJ_con) (L(n_4).Svv_cJ_dag_con)
                       (L(n_4).Cvv_cJ)                 (L(n_4).Rvv_cJ_dag)             (L(n_4).Svv_cJ_con)                     (L(n_4).Pvv_cJ_con - L(n_4).Qvv_cJ_con) ];

L(n_4).Lvv_plusJ =    [(L(n_4).Pvv_pJ + L(n_4).Qvv_pJ) (L(n_4).Svv_pJ)                 (L(n_4).Rvv_pJ)                         (L(n_4).Cvv_pJ_dag)
                       (L(n_4).Svv_pJ_dag)             (L(n_4).Pvv_pJ - L(n_4).Qvv_pJ) (L(n_4).Zvv_pJ)                         (L(n_4).Rvv_pJ)
                       (L(n_4).Rvv_pJ_dag)             (L(n_4).Zvv_pJ_dag)             (L(n_4).Pvv_pJ_con - L(n_4).Qvv_pJ_con) (L(n_4).Svv_pJ_dag_con)
                       (L(n_4).Cvv_pJ)                 (L(n_4).Rvv_pJ_dag)             (L(n_4).Svv_pJ_con)                     (L(n_4).Pvv_pJ_con - L(n_4).Qvv_pJ_con) ];

L(n_4).Lvv_minusJ =   [(L(n_4).Pvv_mJ + L(n_4).Qvv_mJ) (L(n_4).Svv_mJ)                 (L(n_4).Rvv_mJ)                         (L(n_4).Cvv_mJ_dag)
                       (L(n_4).Svv_mJ_dag)             (L(n_4).Pvv_mJ - L(n_4).Qvv_mJ) (L(n_4).Zvv_mJ)                         (L(n_4).Rvv_mJ)
                       (L(n_4).Rvv_mJ_dag)             (L(n_4).Zvv_mJ_dag)             (L(n_4).Pvv_mJ_con - L(n_4).Qvv_mJ_con) (L(n_4).Svv_mJ_dag_con)
                       (L(n_4).Cvv_mJ)                 (L(n_4).Rvv_mJ_dag)             (L(n_4).Svv_mJ_con)                     (L(n_4).Pvv_mJ_con - L(n_4).Qvv_mJ_con) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%%%% Valence- Spin orbit band block (G8p:G7p)
%%%% Interface
L(n_4).Lvs_centerI =  [(sqrt(3/2)*L(n_4).SIGvs_cI_con) (sqrt(2)*L(n_4).Qvs_cI)                 
                       (sqrt(2)*L(n_4).Rvs_cI_con)     (sqrt(1/2)*L(n_4).Svs_cI_con)
                       (sqrt(1/2)*L(n_4).Svs_cI)       (-sqrt(2)*L(n_4).Rvs_cI)         
                       (-sqrt(2)*L(n_4).Qvs_cI)        (sqrt(3/2)*L(n_4).SIGvs_cI)  ];
                   
L(n_4).Lvs_centerIR =  [(sqrt(3/2)*L(n_4).SIGvs_cIR_con) (sqrt(2)*L(n_4).Qvs_cIR)                 
                       (sqrt(2)*L(n_4).Rvs_cIR_con)     (sqrt(1/2)*L(n_4).Svs_cIR_con)
                       (sqrt(1/2)*L(n_4).Svs_cIR)       (-sqrt(2)*L(n_4).Rvs_cIR)         
                       (-sqrt(2)*L(n_4).Qvs_cIR)        (sqrt(3/2)*L(n_4).SIGvs_cIR)  ];
                   
L(n_4).Lvs_centerIL =  [(sqrt(3/2)*L(n_4).SIGvs_cIL_con) (sqrt(2)*L(n_4).Qvs_cIL)                 
                       (sqrt(2)*L(n_4).Rvs_cIL_con)     (sqrt(1/2)*L(n_4).Svs_cIL_con)
                       (sqrt(1/2)*L(n_4).Svs_cIL)       (-sqrt(2)*L(n_4).Rvs_cIL)         
                       (-sqrt(2)*L(n_4).Qvs_cIL)        (sqrt(3/2)*L(n_4).SIGvs_cIL)  ];
                   
L(n_4).Lvs_plusI  =   [(sqrt(3/2)*L(n_4).SIGvs_pI_con) (sqrt(2)*L(n_4).Qvs_pI)                 
                       (sqrt(2)*L(n_4).Rvs_pI_con)     (sqrt(1/2)*L(n_4).Svs_pI_con)
                       (sqrt(1/2)*L(n_4).Svs_pI)       (-sqrt(2)*L(n_4).Rvs_pI)         
                       (-sqrt(2)*L(n_4).Qvs_pI)        (sqrt(3/2)*L(n_4).SIGvs_pI)  ];

L(n_4).Lvs_minusI =  [(sqrt(3/2)*L(n_4).SIGvs_mI_con) (sqrt(2)*L(n_4).Qvs_mI)                 
                       (sqrt(2)*L(n_4).Rvs_mI_con)     (sqrt(1/2)*L(n_4).Svs_mI_con)
                       (sqrt(1/2)*L(n_4).Svs_mI)       (-sqrt(2)*L(n_4).Rvs_mI)         
                       (-sqrt(2)*L(n_4).Qvs_mI)        (sqrt(3/2)*L(n_4).SIGvs_mI)  ];
                   
%%%% LAYER                    
L(n_4).Lvs_centerJ =  [(sqrt(3/2)*L(n_4).SIGvs_cJ_con) (sqrt(2)*L(n_4).Qvs_cJ)                 
                       (sqrt(2)*L(n_4).Rvs_cJ_con)     (sqrt(1/2)*L(n_4).Svs_cJ_con)
                       (sqrt(1/2)*L(n_4).Svs_cJ)       (-sqrt(2)*L(n_4).Rvs_cJ)         
                       (-sqrt(2)*L(n_4).Qvs_cJ)        (sqrt(3/2)*L(n_4).SIGvs_cJ)  ];
                                      
L(n_4).Lvs_plusJ  =   [(sqrt(3/2)*L(n_4).SIGvs_pJ_con) (sqrt(2)*L(n_4).Qvs_pJ)                 
                       (sqrt(2)*L(n_4).Rvs_pJ_con)     (sqrt(1/2)*L(n_4).Svs_pJ_con)
                       (sqrt(1/2)*L(n_4).Svs_pJ)       (-sqrt(2)*L(n_4).Rvs_pJ)         
                       (-sqrt(2)*L(n_4).Qvs_pJ)        (sqrt(3/2)*L(n_4).SIGvs_pJ)  ];

L(n_4).Lvs_minusJ =  [(sqrt(3/2)*L(n_4).SIGvs_mJ_con) (sqrt(2)*L(n_4).Qvs_mJ)                 
                       (sqrt(2)*L(n_4).Rvs_mJ_con)     (sqrt(1/2)*L(n_4).Svs_mJ_con)
                       (sqrt(1/2)*L(n_4).Svs_mJ)       (-sqrt(2)*L(n_4).Rvs_mJ)         
                       (-sqrt(2)*L(n_4).Qvs_mJ)        (sqrt(3/2)*L(n_4).SIGvs_mJ)  ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%%%% Valence- Spin orbit band block (G8p:G7p)*
%%%% Interface
L(n_4).Lsv_centerI  =  [(sqrt(3/2)*L(n_4).SIGvs_cI_dag_con) (sqrt(2)*L(n_4).Rvs_cI_dag_con)   (sqrt(1/2)*L(n_4).Svs_cI_dag) (-sqrt(2)*L(n_4).Qvs_cI_dag)         
                        (sqrt(2)*L(n_4).Qvs_cI_dag)         (sqrt(1/2)*L(n_4).Svs_cI_dag_con) (-sqrt(2)*L(n_4).Rvs_cI_dag)  (sqrt(3/2)*L(n_4).SIGvs_cI_dag)];
                   
L(n_4).Lsv_centerIR =  [(sqrt(3/2)*L(n_4).SIGvs_cIR_dag_con) (sqrt(2)*L(n_4).Rvs_cIR_dag_con)   (sqrt(1/2)*L(n_4).Svs_cIR_dag) (-sqrt(2)*L(n_4).Qvs_cIR_dag)         
                        (sqrt(2)*L(n_4).Qvs_cIR_dag)         (sqrt(1/2)*L(n_4).Svs_cIR_dag_con) (-sqrt(2)*L(n_4).Rvs_cIR_dag)  (sqrt(3/2)*L(n_4).SIGvs_cIR_dag)];
                   
L(n_4).Lsv_centerI  =  [(sqrt(3/2)*L(n_4).SIGvs_cIL_dag_con) (sqrt(2)*L(n_4).Rvs_cIL_dag_con)   (sqrt(1/2)*L(n_4).Svs_cIL_dag) (-sqrt(2)*L(n_4).Qvs_cIL_dag)         
                        (sqrt(2)*L(n_4).Qvs_cIL_dag)         (sqrt(1/2)*L(n_4).Svs_cIL_dag_con) (-sqrt(2)*L(n_4).Rvs_cIL_dag)  (sqrt(3/2)*L(n_4).SIGvs_cIL_dag)];

L(n_4).Lsv_plusI    =  [(sqrt(3/2)*L(n_4).SIGvs_pI_dag_con) (sqrt(2)*L(n_4).Rvs_pI_dag_con)   (sqrt(1/2)*L(n_4).Svs_pI_dag) (-sqrt(2)*L(n_4).Qvs_pI_dag)         
                        (sqrt(2)*L(n_4).Qvs_pI_dag)         (sqrt(1/2)*L(n_4).Svs_pI_dag_con) (-sqrt(2)*L(n_4).Rvs_pI_dag)  (sqrt(3/2)*L(n_4).SIGvs_pI_dag)];

L(n_4).Lsv_minusI   =  [(sqrt(3/2)*L(n_4).SIGvs_mI_dag_con) (sqrt(2)*L(n_4).Rvs_mI_dag_con)   (sqrt(1/2)*L(n_4).Svs_mI_dag) (-sqrt(2)*L(n_4).Qvs_mI_dag)         
                        (sqrt(2)*L(n_4).Qvs_mI_dag)         (sqrt(1/2)*L(n_4).Svs_mI_dag_con) (-sqrt(2)*L(n_4).Rvs_mI_dag)  (sqrt(3/2)*L(n_4).SIGvs_mI_dag)];
                   
%%%% LAYER                    
L(n_4).Lsv_centerJ  =  [(sqrt(3/2)*L(n_4).SIGvs_cJ_dag_con) (sqrt(2)*L(n_4).Rvs_cJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_cJ_dag) (-sqrt(2)*L(n_4).Qvs_cJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_cJ_dag)         (sqrt(1/2)*L(n_4).Svs_cJ_dag_con) (-sqrt(2)*L(n_4).Rvs_cJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_cJ_dag)];

L(n_4).Lsv_plusJ    =  [(sqrt(3/2)*L(n_4).SIGvs_pJ_dag_con) (sqrt(2)*L(n_4).Rvs_pJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_pJ_dag) (-sqrt(2)*L(n_4).Qvs_pJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_pJ_dag)         (sqrt(1/2)*L(n_4).Svs_pJ_dag_con) (-sqrt(2)*L(n_4).Rvs_pJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_pJ_dag)];

L(n_4).Lsv_minusJ   =  [(sqrt(3/2)*L(n_4).SIGvs_mJ_dag_con) (sqrt(2)*L(n_4).Rvs_mJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_mJ_dag) (-sqrt(2)*L(n_4).Qvs_mJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_mJ_dag)         (sqrt(1/2)*L(n_4).Svs_mJ_dag_con) (-sqrt(2)*L(n_4).Rvs_mJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_mJ_dag)];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Spin orbit band block
%%%% Interface
L(n_4).Lss_centerI =  [(L(n_4).Ps_cI)     (L(n_4).Cs_cI) 
                       (L(n_4).Cs_cI_dag) (L(n_4).Ps_cI_con) ];

L(n_4).Lss_centerIR =  [(L(n_4).Ps_cIR)     (L(n_4).Cs_cIR) 
                       (L(n_4).Cs_cIR_dag)  (L(n_4).Ps_cIR_con) ];                   
                   
L(n_4).Lss_centerIL =  [(L(n_4).Ps_cIL)     (L(n_4).Cs_cIL) 
                       (L(n_4).Cs_cIL_dag)  (L(n_4).Ps_cIL_con) ];
                   
L(n_4).Lss_plusI   =  [(L(n_4).Ps_pI)     (L(n_4).Cs_pI) 
                       (L(n_4).Cs_pI_dag) (L(n_4).Ps_pI_con) ];

L(n_4).Lss_minusI  =  [(L(n_4).Ps_mI)     (L(n_4).Cs_mI) 
                       (L(n_4).Cs_mI_dag) (L(n_4).Ps_mI_con) ];
%%%% LAYER                    
L(n_4).Lss_centerJ =  [(L(n_4).Ps_cJ)     (L(n_4).Cs_cJ) 
                       (L(n_4).Cs_cJ_dag) (L(n_4).Ps_cJ_con) ];

L(n_4).Lss_plusJ   =  [(L(n_4).Ps_pJ)     (L(n_4).Cs_pJ) 
                       (L(n_4).Cs_pJ_dag) (L(n_4).Ps_pJ_con) ];

L(n_4).Lss_minusJ  =  [(L(n_4).Ps_mJ)     (L(n_4).Cs_mJ) 
                       (L(n_4).Cs_mJ_dag) (L(n_4).Ps_mJ_con) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Construct Hamiltonian that will be solved using finite difference
%%%% INTERFACE 

if Bands == 6

L(n_4).H_centerI = [  ((L(n_4).Lvv_centerI)+(L(n_4).Evv_centerI)) (L(n_4).Lvs_centerI)
                      (L(n_4).Lsv_centerI) ((L(n_4).Lss_centerI)+(L(n_4).Ess_centerI)) ];

L(n_4).H_centerIR= [  ((L(n_4).Lvv_centerIR)+(L(n_4).Evv_centerIR)) (L(n_4).Lvs_centerIR)
                      (L(n_4).Lsv_centerIR) (L(n_4).Lss_centerIR)+(L(n_4).Ess_centerIR)) ];

L(n_4).H_centerIL= [  ((L(n_4).Lvv_centerIL)+(L(n_4).Evv_centerIL)) (L(n_4).Lvs_centerIL)
                      (L(n_4).Lsv_centerIL) (L(n_4).Lss_centerIL)+(L(n_4).Ess_centerIL)) ];

L(n_4).H_plusI   = [ (L(n_4).Lvv_plusI) (L(n_4).Lvs_plusI)
                     (L(n_4).Lsv_plusI) (L(n_4).Lss_plusI) ];
                 
L(n_4).H_minusI  = [ (L(n_4).Lvv_minusI) (L(n_4).Lvs_minusI)
                     (L(n_4).Lsv_minusI) (L(n_4).Lss_minusI) ];
%%%% LAYER 
L(n_4).H_centerJ = [ ((L(n_4).Lvv_centerJ)+(L(n_4).Evv_centerJ)) (L(n_4).Lvs_centerJ)
                     (L(n_4).Lsv_centerJ) ((L(n_4).Lss_centerJ)+(L(n_4).Ess_centerJ)) ];

L(n_4).H_plusJ   = [ (L(n_4).Lvv_plusJ) (L(n_4).Lvs_plusJ)
                     (L(n_4).Lsv_plusJ) (L(n_4).Lss_plusJ) ];
                 
L(n_4).H_minusJ  = [ (L(n_4).Lvv_minusJ) (L(n_4).Lvs_minusJ)
                     (L(n_4).Lsv_minusJ) (L(n_4).Lss_minusJ) ];
elseif Bands == 8
% L(n_4).H_centerI = [ (L(n_4).Lcc_centerI) (L(n_4).Kcv_centerI) (zeros(2,2))
%                      (L(n_4).Kvc_centerI) (L(n_4).Lvv_centerI) (L(n_4).Lvs_centerI)
%                      (zeros(2,2))         (L(n_4).Lsv_centerI) (L(n_4).Lss_centerI) ];
% 
% L(n_4).H_centerIR= [ (L(n_4).Lcc_centerIR) (L(n_4).Kcv_centerIR) (zeros(2,2))
%                      (L(n_4).Kvc_centerIR) (L(n_4).Lvv_centerIR) (L(n_4).Lvs_centerIR)
%                      (zeros(2,2))          (L(n_4).Lsv_centerIR) (L(n_4).Lss_centerIR) ];
% 
% L(n_4).H_centerIL= [ (L(n_4).Lcc_centerIL) (L(n_4).Kcv_centerIL) (zeros(2,2))
%                      (L(n_4).Kvc_centerIL) (L(n_4).Lvv_centerIL) (L(n_4).Lvs_centerIL)
%                      (zeros(2,2))          (L(n_4).Lsv_centerIL) (L(n_4).Lss_centerIL) ];
% 
% L(n_4).H_plusI   = [ (L(n_4).Lcc_plusI) (L(n_4).Kcv_plusI) (zeros(2,2))
%                      (L(n_4).Kvc_plusI) (L(n_4).Lvv_plusI) (L(n_4).Lvs_plusI)
%                      (zeros(2,2))       (L(n_4).Lsv_plusI) (L(n_4).Lss_plusI) ];
%                  
% L(n_4).H_minusI  = [ (L(n_4).Lcc_minusI) (L(n_4).Kcv_minusI) (zeros(2,2))
%                      (L(n_4).Kvc_minusI) (L(n_4).Lvv_minusI) (L(n_4).Lvs_minusI)
%                      (zeros(2,2))        (L(n_4).Lsv_minusI) (L(n_4).Lss_minusI) ];
% %%%% LAYER 
% L(n_4).H_centerJ = [ (L(n_4).Lcc_centerJ) (L(n_4).Kcv_centerJ) (zeros(2,2))
%                      (L(n_4).Kvc_centerJ) (L(n_4).Lvv_centerJ) (L(n_4).Lvs_centerJ)
%                      (zeros(2,2))         (L(n_4).Lsv_centerJ) (L(n_4).Lss_centerJ) ];
% 
% L(n_4).H_plusJ   = [ (L(n_4).Lcc_plusJ) (L(n_4).Kcv_plusJ) (zeros(2,2))
%                      (L(n_4).Kvc_plusJ) (L(n_4).Lvv_plusJ) (L(n_4).Lvs_plusJ)
%                      (zeros(2,2))       (L(n_4).Lsv_plusJ) (L(n_4).Lss_plusJ) ];
%                  
% L(n_4).H_minusJ  = [ (L(n_4).Lcc_minusJ) (L(n_4).Kcv_minusJ) (zeros(2,2))
%                      (L(n_4).Kvc_minusJ) (L(n_4).Lvv_minusJ) (L(n_4).Lvs_minusJ)
%                      (zeros(2,2))        (L(n_4).Lsv_minusJ) (L(n_4).Lss_minusJ) ];    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%% INTERFACE   
%     L(n_4).H11centerI =  (L(n_4).E_elI + (L(n_4).M_ecI));
%     L(n_4).H11plusI   =  (L(n_4).M_epI);         
%     L(n_4).H11minusI  =  (L(n_4).M_emI);
% %%%% LAYER   
%     L(n_4).H11centerJ =  (L(n_4).E_elJ + (L(n_4).M_ecJ));
%     L(n_4).H11plusJ   =  (L(n_4).M_epJ);         
%     L(n_4).H11minusJ  =  (L(n_4).M_emJ);
%     
% %Define the Hamiltonian that will be solved in the finite difference
% %theorum Holes
% %%%% INTERFACE POINT LEFT
%     L(n_4).H33centerIL =  [(L(n_4).E_hhIL+(L(n_4).PcIL+(L(n_4).QcIL-L(n_4).Qe_shearIL)))  (0)                                                       (0)                                                        (L(n_4).RcIL)                                               (0)                                          (sqrt(2)*L(n_4).RcIL)
%                           (0)                                                        (L(n_4).E_hhIL+(L(n_4).PcIL+(L(n_4).QcIL-L(n_4).Qe_shearIL))) (-L(n_4).RcIL_dag)                                          (0)                                                        (-sqrt(2)*L(n_4).RcIL_dag)                    (0)
%                           (0)                                                        (-L(n_4).RcIL)                                             (L(n_4).E_lhIL+(L(n_4).PcIL-(L(n_4).QcIL-L(n_4).Qe_shearIL)))  (0)                                                        (sqrt(2)*(L(n_4).QcIL - L(n_4).Qe_shearIL))    (0)
%                           (L(n_4).RcIL_dag)                                           (0)                                                       (0)                                                        (L(n_4).E_lhIL+(L(n_4).PcIL-(L(n_4).QcIL-L(n_4).Qe_shearIL)))  (0)                                          (sqrt(2)*(L(n_4).QcIL - L(n_4).Qe_shearIL)) 
%                           (0)                                                        (-sqrt(2)*L(n_4).RcIL)                                     (sqrt(2)*(L(n_4).QcIL -L(n_4).Qe_shearIL))                   (0)                                                        (L(n_4).E_soIL + L(n_4).PcIL - L(n_4).DeltaIL)  (0)
%                           (sqrt(2)*L(n_4).RcIL_dag)                                   (0)                                                       (0)                                                        (sqrt(2)*(L(n_4).QcIL -L(n_4).Qe_shearIL))                   (0)                                          (L(n_4).E_soIL + L(n_4).PcIL - L(n_4).DeltaIL)    ];
% 
% %%%% INTERFACE POINT RIGHT
%     L(n_4).H33centerIR =  [(L(n_4).E_hhIR+(L(n_4).PcIR+(L(n_4).QcIR-L(n_4).Qe_shearIR)))  (0)                                                       (0)                                                        (L(n_4).RcIR)                                               (0)                                          (sqrt(2)*L(n_4).RcIR)
%                           (0)                                                        (L(n_4).E_hhIR+(L(n_4).PcIR+(L(n_4).QcIR-L(n_4).Qe_shearIR))) (-L(n_4).RcIR_dag)                                          (0)                                                        (-sqrt(2)*L(n_4).RcIR_dag)                    (0)
%                           (0)                                                        (-L(n_4).RcIR)                                             (L(n_4).E_lhIR+(L(n_4).PcIR-(L(n_4).QcIR-L(n_4).Qe_shearIR)))  (0)                                                        (sqrt(2)*(L(n_4).QcIR - L(n_4).Qe_shearIR))    (0)
%                           (L(n_4).RcIR_dag)                                           (0)                                                       (0)                                                        (L(n_4).E_lhIR+(L(n_4).PcIR-(L(n_4).QcIR-L(n_4).Qe_shearIR)))  (0)                                          (sqrt(2)*(L(n_4).QcIR - L(n_4).Qe_shearIR)) 
%                           (0)                                                        (-sqrt(2)*L(n_4).RcIR)                                     (sqrt(2)*(L(n_4).QcIR -L(n_4).Qe_shearIR))                   (0)                                                        (L(n_4).E_soIR + L(n_4).PcIR - L(n_4).DeltaIR)  (0)
%                           (sqrt(2)*L(n_4).RcIR_dag)                                   (0)                                                       (0)                                                        (sqrt(2)*(L(n_4).QcIR -L(n_4).Qe_shearIR))                   (0)                                          (L(n_4).E_soIR + L(n_4).PcIR - L(n_4).DeltaIR)    ];
% 
% 
% 
% %%%% INTERFACE POINT CENTER
%     L(n_4).H33centerI =  [(L(n_4).E_hhI+(L(n_4).PcI+(L(n_4).QcI-L(n_4).Qe_shearI)))  (0)                                                       (0)                                                        (L(n_4).RcI)                                               (0)                                          (sqrt(2)*L(n_4).RcI)
%                           (0)                                                        (L(n_4).E_hhI+(L(n_4).PcI+(L(n_4).QcI-L(n_4).Qe_shearI))) (-L(n_4).RcI_dag)                                          (0)                                                        (-sqrt(2)*L(n_4).RcI_dag)                    (0)
%                           (0)                                                        (-L(n_4).RcI)                                             (L(n_4).E_lhI+(L(n_4).PcI-(L(n_4).QcI-L(n_4).Qe_shearI)))  (0)                                                        (sqrt(2)*(L(n_4).QcI - L(n_4).Qe_shearI))    (0)
%                           (L(n_4).RcI_dag)                                           (0)                                                       (0)                                                        (L(n_4).E_lhI+(L(n_4).PcI-(L(n_4).QcI-L(n_4).Qe_shearI)))  (0)                                          (sqrt(2)*(L(n_4).QcI - L(n_4).Qe_shearI)) 
%                           (0)                                                        (-sqrt(2)*L(n_4).RcI)                                     (sqrt(2)*(L(n_4).QcI -L(n_4).Qe_shearI))                   (0)                                                        (L(n_4).E_soI + L(n_4).PcI - L(n_4).DeltaI)  (0)
%                           (sqrt(2)*L(n_4).RcI_dag)                                   (0)                                                       (0)                                                        (sqrt(2)*(L(n_4).QcI -L(n_4).Qe_shearI))                   (0)                                          (L(n_4).E_soI + L(n_4).PcI - L(n_4).DeltaI)    ];
% 
%     L(n_4).H33plusI   =  [(L(n_4).PpI + L(n_4).QpI)             (0)                                (-1i)*-(L(n_4).SpI*k_m)                    (0)                                         (-1i)*((L(n_4).SpI)*k_m/sqrt(2))        (0)
%                           (0)                                   (L(n_4).PpI + L(n_4).QpI)          (0)                                        (-1i)*(-L(n_4).SpI*k_p)                     (0)                                     (-1i)*((L(n_4).SpI)*k_p/sqrt(2))
%                           (-1i)*((-L(n_4).SpI_dag)*k_p)         (0)                                (L(n_4).PpI - L(n_4).QpI)                  (-1i)*(L(n_4).CpI*k_m)                      (sqrt(2)*L(n_4).QpI)                    (-1i)*(sqrt(3/2)*L(n_4).EpsilonpI*k_m)
%                           (0)                                   (-1i)*((-L(n_4).SpI_dag)*k_m)      (-1i)*(L(n_4).CpI_dag*k_p)                 (L(n_4).PpI - L(n_4).QpI)                   (-1i)*(-sqrt(3/2)*L(n_4).EpsilonpI*k_p) (sqrt(2)*L(n_4).QpI)
%                           (-1i)*((L(n_4).SpI_dag)*k_p/sqrt(2))  (0)                                (sqrt(2)*L(n_4).QpI)                       (-1i)*(-sqrt(3/2)*L(n_4).EpsilonpI_dag*k_m) (L(n_4).PpI)                            (-1i)*(-L(n_4).CpI*k_m)                      
%                           (0)                                   (-1i)*(L(n_4).SpI_dag*k_m/sqrt(2)) (-1i)*(sqrt(3/2)*L(n_4).EpsilonpI_dag*k_p) (sqrt(2)*L(n_4).QpI)                        (-1i)*(-L(n_4).CpI_dag*k_p)             (L(n_4).PpI)                 ];
%                
%     L(n_4).H33minusI  =  [(L(n_4).PmI + L(n_4).QmI)          (0)                                (-1i)*(-L(n_4).SmI*k_m)                    (0)                                         (-1i)*((L(n_4).SmI)*k_m/sqrt(2))          (0)
%                           (0)                                (L(n_4).PmI + L(n_4).QmI)          (0)                                        (-1i)*(-L(n_4).SmI*k_p)                     (0)                                     (-1i)*((L(n_4).SmI)*k_p/sqrt(2))
%                           (-1i)*((-L(n_4).SmI_dag)*k_p)        (0)                                (L(n_4).PmI - L(n_4).QmI)                  (-1i)*(L(n_4).CmI*k_m)                      (sqrt(2)*L(n_4).QmI)                    (-1i)*(sqrt(3/2)*L(n_4).EpsilonmI*k_m)
%                           (0)                                (-1i)*((-L(n_4).SmI_dag)*k_m)        (-1i)*(L(n_4).CmI_dag*k_p)                 (L(n_4).PmI - L(n_4).QmI)                   (-1i)*(-sqrt(3/2)*L(n_4).EpsilonmI*k_p) (sqrt(2)*L(n_4).QmI)
%                           (-1i)*((L(n_4).SmI_dag)*k_p/sqrt(2)) (0)                                (sqrt(2)*L(n_4).QmI)                       (-1i)*(-sqrt(3/2)*L(n_4).EpsilonmI_dag*k_m) (L(n_4).PmI)                            (-1i)*(-L(n_4).CmI*k_m)                      
%                           (0)                                (-1i)*((L(n_4).SmI_dag)*k_m/sqrt(2)) (-1i)*(sqrt(3/2)*L(n_4).EpsilonmI_dag*k_p) (sqrt(2)*L(n_4).QmI)                        (-1i)*(-L(n_4).CmI_dag*k_p)             (L(n_4).PmI)                 ];                
%                       
%  %%%% LAYER   
%     L(n_4).H33centerJ =  [(L(n_4).E_hhJ+(L(n_4).PcJ+(L(n_4).QcJ-L(n_4).Qe_shearJ)))  (0)                                                       (0)                                                        (L(n_4).RcJ)                                               (0)                                          (sqrt(2)*L(n_4).RcJ)
%                           (0)                                                        (L(n_4).E_hhJ+(L(n_4).PcJ+(L(n_4).QcJ-L(n_4).Qe_shearJ))) (-L(n_4).RcJ_dag)                                          (0)                                                        (-sqrt(2)*L(n_4).RcJ_dag)                    (0)
%                           (0)                                                        (-L(n_4).RcJ)                                             (L(n_4).E_lhJ+(L(n_4).PcJ-(L(n_4).QcJ-L(n_4).Qe_shearJ)))  (0)                                                        (sqrt(2)*(L(n_4).QcJ - L(n_4).Qe_shearJ))    (0)
%                           (L(n_4).RcJ_dag)                                           (0)                                                       (0)                                                        (L(n_4).E_lhJ+(L(n_4).PcJ-(L(n_4).QcJ-L(n_4).Qe_shearJ)))  (0)                                          (sqrt(2)*(L(n_4).QcJ - L(n_4).Qe_shearJ)) 
%                           (0)                                                        (-sqrt(2)*L(n_4).RcJ)                                     (sqrt(2)*(L(n_4).QcJ -L(n_4).Qe_shearJ))                   (0)                                                        (L(n_4).E_soJ + L(n_4).PcJ - L(n_4).DeltaJ)  (0)
%                           (sqrt(2)*L(n_4).RcJ_dag)                                   (0)                                                       (0)                                                        (sqrt(2)*(L(n_4).QcJ -L(n_4).Qe_shearJ))                   (0)                                          (L(n_4).E_soJ + L(n_4).PcJ - L(n_4).DeltaJ)    ];
% 
%     L(n_4).H33plusJ   =  [(L(n_4).PpJ + L(n_4).QpJ)    	 (0)                            (-1i)*-(L(n_4).SpJ*k_m)          (0)                               (-1i)*(L(n_4).SpJ*k_m/sqrt(2))    (0)
%                           (0)                          	 (L(n_4).PpJ + L(n_4).QpJ)      (0)                              (-1i)*(-L(n_4).SpJ*k_p)           (0)                               (-1i)*(L(n_4).SpJ*k_p/sqrt(2))
%                           (-1i)*-(L(n_4).SpJ*k_p)        (0)                            (L(n_4).PpJ - L(n_4).QpJ)        (0)                               (sqrt(2)*L(n_4).QpJ)              (-1i)*(sqrt(3/2)*L(n_4).EpsilonpJ*k_m)
%                           (0)                            (-1i)*-(L(n_4).SpJ*k_m)        (0)                              (L(n_4).PpJ - L(n_4).QpJ)         (-1i)*(-sqrt(3/2)*L(n_4).EpsilonpJ*k_p) (sqrt(2)*L(n_4).QpJ)
%                           (-1i)*(L(n_4).SpJ*k_p/sqrt(2)) (0)                            (sqrt(2)*L(n_4).QpJ)             (-1i)*(-sqrt(3/2)*L(n_4).EpsilonpJ*k_m) (L(n_4).PpJ)                      (0)                      
%                           (0)                            (-1i)*(L(n_4).SpJ*k_m/sqrt(2)) (-1i)*(sqrt(3/2)*L(n_4).EpsilonpJ*k_p) (sqrt(2)*L(n_4).QpJ)              (0)                               (L(n_4).PpJ)                 ];
%                  
%     L(n_4).H33minusJ  =  [(L(n_4).PmJ + L(n_4).QmJ)      (0)                            (-1i)*(-L(n_4).SmJ*k_m)          (0)                               (-1i)*(L(n_4).SmJ*k_m/sqrt(2))    (0)
%                           (0)                            (L(n_4).PmJ + L(n_4).QmJ)      (0)                              (-1i)*(-L(n_4).SmJ*k_p)           (0)                               (-1i)*(L(n_4).SmJ*k_p/sqrt(2))
%                           (-1i)*(-L(n_4).SmJ*k_p)        (0)                            (L(n_4).PmJ - L(n_4).QmJ)        (0)                    	       (sqrt(2)*L(n_4).QmJ)              (-1i)*(sqrt(3/2)*L(n_4).EpsilonmJ*k_m)
%                           (0)                            (-1i)*(-L(n_4).SmJ*k_m)        (0)                              (L(n_4).PmJ - L(n_4).QmJ)         (-1i)*(-sqrt(3/2)*L(n_4).EpsilonmJ*k_p) (sqrt(2)*L(n_4).QmJ)
%                           (-1i)*(L(n_4).SmJ*k_p/sqrt(2)) (0)                            (sqrt(2)*L(n_4).QmJ)             (-1i)*(-sqrt(3/2)*L(n_4).EpsilonmJ*k_m) (L(n_4).PmJ)                      (0)                      
%                           (0)                            (-1i)*(L(n_4).SmJ*k_m/sqrt(2)) (-1i)*(sqrt(3/2)*L(n_4).EpsilonmJ*k_p) (sqrt(2)*L(n_4).QmJ)              (0)                               (L(n_4).PmJ)                 ];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for n_5 = 1 : n_layer
A=L(n_5).layer_npts;
%%%% Electrons
% Main_diagonal1_e = repmat(L(n_5).H11centerJ,[1,A-1]);
% Main_diagonal2_e = [ (L(n_5).H11centerI) , Main_diagonal1_e];
% Super_diagonal1_e = repmat(L(n_5).H11plusJ,[1,A-1]);
% Super_diagonal2_e = [ (L(n_5).H11plusI) , Super_diagonal1_e];
% Sub_diagonal1_e = repmat(L(n_5).H11minusJ,[1,A-1]); 
% Sub_diagonal2_e = [ (L(n_5).H11minusI) , Sub_diagonal1_e];
% 
% Main_diag_e{n_5} = Main_diagonal2_e;
% Super_diag_e{n_5}= Super_diagonal2_e ;
% Sub_diag_e{n_5}  = Sub_diagonal2_e ;
%%%% Holes

if (n_5 == 1)
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2_h = [ (L(n_5).H_centerI) , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-1]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-1]); 
Sub_diagonal2_h = [ (L(n_5).H_minusI) , Sub_diagonal1_h];

elseif (n_5 == n_layer)
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2_h = [ (L(n_5).H_centerI) , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-1]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-1]); 
Sub_diagonal2_h = [ (L(n_5).H_minusI) , Sub_diagonal1_h];

elseif (L(n_5+1).layer_cntr == 1) 
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-2]);
Main_diagonal2_h = [ (L(n_5).H_centerI) , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-2]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h , (L(n_5+1).H33minusI)' ];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-1]); 
Sub_diagonal2_h = [ (L(n_5).H_minusI) , Sub_diagonal1_h];

elseif (L(n_5).layer_cntr == 1)
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-3]);
Main_diagonal2_h = [L(n_5).H_centerIL , (L(n_5).H33centerI) , L(n_5).H33centerIR , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-2]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h , (L(n_5+1).H33minusI)' ];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-2]); 
Sub_diagonal2_h = [ L(n_5).H_minusI , (L(n_5).H33plusI)' , Sub_diagonal1_h ];

elseif (L(n_5-1).layer_cntr == 1)
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-2]);
Main_diagonal2_h = [L(n_5).H_centerIL , (L(n_5).H33centerI) , L(n_5).H33centerIR , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-1]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-2]); 
Sub_diagonal2_h = [ (L(n_5).H_minusI) , (L(n_5).H33plusI)' , Sub_diagonal1_h]; 

else
Main_diagonal1_h = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2_h = [ (L(n_5).H_centerI) , Main_diagonal1_h];
Super_diagonal1_h = repmat(L(n_5).H_plusJ,[1,A-1]);
Super_diagonal2_h = [ (L(n_5).H_plusI) , Super_diagonal1_h];
Sub_diagonal1_h = repmat(L(n_5).H_minusJ,[1,A-1]); 
Sub_diagonal2_h = [ (L(n_5).H_minusI) , Sub_diagonal1_h];
end   
Main_diag_h{n_5} = Main_diagonal2_h;
Super_diag_h{n_5}= Super_diagonal2_h ;
Sub_diag_h{n_5}  = Sub_diagonal2_h ;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert the N cell array consisting of diagonal, sub and super 
%%%% diagonal terms into matrices
%%%% Electrons
% Main_diag_e1=cell2mat(Main_diag_e);
% Main_diag_e2= [Main_diag_e1 , infboundary_e1];
% Main_diag_e3= V_EfieldE + Main_diag_e2;
% Super_diag_e1=cell2mat(Super_diag_e);
% Super_diag_e1(:,tot_mesh_pts)=[];
% Super_diag_e2= [ Super_diag_e1 , infboundary_e2];
% Sub_diag_e1=cell2mat(Sub_diag_e);
% Sub_diag_e1(:,1)=[];
% Sub_diag_e1(:,1)=0;
% Sub_diag_e2= [Sub_diag_e1 , infboundary_e2];
% Hmatrix_Electrons = sparse(diag(Main_diag_e3,0) + diag(Super_diag_e2,1) + diag(Sub_diag_e2,-1)); 
%%%% Holes
%%%% Re-shape the N by 6 matrices into N, 6 by 6 matrices.
Amd1 =cell2mat(Main_diag_h);
Amd2 = [Amd1  , infboundary_h1];
Asup1 = cell2mat(Super_diag_h);
Asup2 = [Asup1  , infboundary_h2];
Asup2(:,((6*tot_mesh_pts-5):6*tot_mesh_pts))=[];
Asub1 = cell2mat(Sub_diag_h);
Asub1(:,7:12)=[];
Asub2 = [Asub1  , infboundary_h2];

Amd = reshape(Amd2,[6,6,tot_mesh_pts+1]);
Asup = reshape(Asup2,[6,6,(tot_mesh_pts)]);
Asub = reshape(Asub2,[6,6,(tot_mesh_pts)]);

%%%% Call blktridiag.m which puts the diagonal (Amd), subb-diagonal (Asub)
%%%% and super-diagonal (Asup) 3 by 3 matrices into a block tridiagonal
%%%% matrix B 
Hmatrix_Holes = blktridiag(Amd,Asub,Asup) + sparse(V_EfieldH); 
setup_time=setup_time+toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.5 Use the function eigs to calculate eigenvalues and 
%%%% eigenvectors for Conduction and Valence bands. 
%%%% ELECTRONS
%[V_valelectrons,E_valelectrons] = eigs(Hmatrix_Electrons,n_Esub, 0.1, options);
%%%% HOLES
[V_valholes,E_valholes, flag] = eigs(Hmatrix_Holes ,n_Esub, -0.001, options);
%%%% Sort Eigenvalues and Eigenfunctions
%[Data(Cn_kpts).E_e,Data((Cn_kpts)).index_e] = sort(diag(E_valelectrons),'ascend');
%V_ee = V_valelectrons(:,Data((Cn_kpts)).index_e);
% disp(flag)
[Data(Cn_kpts).E_h,Data((Cn_kpts)).index]  = sort((diag(E_valholes)),'descend');
V_valholes=V_valholes(:,Data((Cn_kpts)).index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if flag_1 == 0
% %%%% Now orders the full envelope functions
% 
% 
% for n_6=1:n_Esub;
% 
%     if mod(n_6,2) == 1
%     Data(Cn_kpts).V_valholesF(:,n_6)=V_valholes(:,n_6);
%     else
%     Data(Cn_kpts).V_valholesF(:,n_6)=(V_valholes(:,n_6)-((V_valholes(:,(n_6-1))'*V_valholes(:,(n_6))).*V_valholes(:,(n_6-1))))/...
%                       sqrt(((V_valholes(:,n_6)-((V_valholes(:,(n_6-1))'*V_valholes(:,(n_6))).*V_valholes(:,(n_6-1))))'...
%                       *(V_valholes(:,n_6)-((V_valholes(:,(n_6-1))'*V_valholes(:,(n_6))).*V_valholes(:,(n_6-1))))));
%     end
% end
% for n_7=1:n_Esub
% vectors=reshape(Data(Cn_kpts).V_valholesF(:,n_7),6,tot_mesh_pts+1)';
% 
% Data(Cn_kpts).V_hhu(:,n_7)=vectors(:,1);
% Data(Cn_kpts).V_hhd(:,n_7)=vectors(:,2);
% Data(Cn_kpts).V_lhu(:,n_7)=vectors(:,3);
% Data(Cn_kpts).V_lhd(:,n_7)=vectors(:,4);
% Data(Cn_kpts).V_sou(:,n_7)=vectors(:,5);
% Data(Cn_kpts).V_sod(:,n_7)=vectors(:,6);
% end
% 
% for CC1=1:n_Esub
% Data(Cn_kpts).Chhu(:,CC1)=Data(Cn_kpts).V_hhu(:,CC1)'*Data(Cn_kpts).V_hhu(:,CC1);
% Data(Cn_kpts).Chhd(:,CC1)=Data(Cn_kpts).V_hhd(:,CC1)'*Data(Cn_kpts).V_hhd(:,CC1);
% Data(Cn_kpts).Clhu(:,CC1)=Data(Cn_kpts).V_lhu(:,CC1)'*Data(Cn_kpts).V_lhu(:,CC1);
% Data(Cn_kpts).Clhd(:,CC1)=Data(Cn_kpts).V_lhd(:,CC1)'*Data(Cn_kpts).V_lhd(:,CC1);
% Data(Cn_kpts).Csou(:,CC1)=Data(Cn_kpts).V_sou(:,CC1)'*Data(Cn_kpts).V_sou(:,CC1);
% Data(Cn_kpts).Csod(:,CC1)=Data(Cn_kpts).V_sod(:,CC1)'*Data(Cn_kpts).V_sod(:,CC1);
% 
% Data(Cn_kpts).C_tot(:,CC1)=Data(Cn_kpts).V_hhu(:,CC1)'*Data(Cn_kpts).V_hhu(:,CC1)+...
%                            Data(Cn_kpts).V_hhd(:,CC1)'*Data(Cn_kpts).V_hhd(:,CC1)+...
%                            Data(Cn_kpts).V_lhu(:,CC1)'*Data(Cn_kpts).V_lhu(:,CC1)+...
%                            Data(Cn_kpts).V_lhd(:,CC1)'*Data(Cn_kpts).V_lhd(:,CC1)+...
%                            Data(Cn_kpts).V_sou(:,CC1)'*Data(Cn_kpts).V_sou(:,CC1)+...
%                            Data(Cn_kpts).V_sod(:,CC1)'*Data(Cn_kpts).V_sod(:,CC1);
% 
% Data(Cn_kpts).Whhu(:,CC1)=abs(Data(Cn_kpts).V_hhu(:,CC1));
% Data(Cn_kpts).Whhd(:,CC1)=abs(Data(Cn_kpts).V_hhd(:,CC1));
% Data(Cn_kpts).Wlhu(:,CC1)=abs(Data(Cn_kpts).V_lhu(:,CC1));
% Data(Cn_kpts).Wlhd(:,CC1)=abs(Data(Cn_kpts).V_lhd(:,CC1));
% Data(Cn_kpts).Wsou(:,CC1)=abs(Data(Cn_kpts).V_sou(:,CC1));
% Data(Cn_kpts).Wsod(:,CC1)=abs(Data(Cn_kpts).V_sod(:,CC1));
% end
% 
% % if Data(Cn_kpts).V_hhu(:,1)'*Data(Cn_kpts).V_hhu(:,1)<Data(Cn_kpts).V_hhu(:,2)'*Data(Cn_kpts).V_hhu(:,2)
% %     Data(Cn_kpts).V_hhu(:,:)=Data(Cn_kpts).V_hhu(:,[2,1]);
% %     Data(Cn_kpts).V_hhd(:,:)=Data(Cn_kpts).V_hhd(:,[2,1]);
% %     Data(Cn_kpts).V_lhu(:,:)=Data(Cn_kpts).V_lhu(:,[2,1]);
% %     Data(Cn_kpts).V_lhd(:,:)=Data(Cn_kpts).V_lhd(:,[2,1]);
% %     Data(Cn_kpts).V_sou(:,:)=Data(Cn_kpts).V_sou(:,[2,1]);
% %     Data(Cn_kpts).V_sod(:,:)=Data(Cn_kpts).V_sod(:,[2,1]);
% %     Data(Cn_kpts).E_h=Data(Cn_kpts).E_h([2;1]);
% % end
% %%%% Renormalise envelope functions at k=0 due to degenerate eigenvalues
% if Cn_kpts == 1
%      Re_norm_u =  Data(Cn_kpts).V_hhu(:,1)'*Data(Cn_kpts).V_hhu(:,1);
%      Data(Cn_kpts).V_hhu(:,1)=(1/sqrt(Re_norm_u)).*(Data(Cn_kpts).V_hhu(:,1));
%      Data(Cn_kpts).V_hhu(:,2)=(0).*(Data(Cn_kpts).V_hhu(:,2));
%      Re_norm_d =  Data(Cn_kpts).V_hhd(:,2)'*Data(Cn_kpts).V_hhd(:,2);
%      Data(Cn_kpts).V_hhd(:,2)=(1/sqrt(Re_norm_d)).*(Data(Cn_kpts).V_hhd(:,2));
%      Data(Cn_kpts).V_hhd(:,1)=(0).*(Data(Cn_kpts).V_hhd(:,1));
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Now calulates the inner products to calculate optical matrix elements
% else
% for n_8=1:n_Esub
% for n_9=1:n_Esub
% vectors=reshape(V_valholes(:,n_9),6,tot_mesh_pts+1)';
% V_hhu(:,n_9)=vectors(:,1);
% V_hhd(:,n_9)=vectors(:,2);
% V_lhu(:,n_9)=vectors(:,3);
% V_lhd(:,n_9)=vectors(:,4);
% V_sou(:,n_9)=vectors(:,5);
% V_sod(:,n_9)=vectors(:,6);
% 
% % %%%% Renormalise envelope functions at k=0 due to degenerate eigenvalues
% % if Cn_kpts == 1
% %      Re_norm_u =  V_hhu(:,1)'*V_hhu(:,1);
% %      V_hhu(:,1)= (1/sqrt(Re_norm_u)).*V_hhu(:,1);
% %      V_hhu(:,2)= (0).*(V_hhu(:,2));
% %      Re_norm_d =  V_hhd(:,2)'*V_hhd(:,2);
% %      V_hhd(:,2)= (1/sqrt(Re_norm_d)).*(V_hhd(:,2));
% %      V_hhd(:,1)= (0).*(V_hhd(:,1));
% % end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Energy seperation E_cv=(Ej-Ev)
% Data(Cn_kpts).E_cv(n_8,n_9)=(Data(Cn_kpts).E_e(n_8)-Data(Cn_kpts).E_h(n_9));
% %%%% Coefficient of <Phi_vi|Phi_vi> for i=HH, LH or SO
% Data(Cn_kpts).chhu(n_9)=V_hhu(:,n_9)'*V_hhu(:,n_9);
% Data(Cn_kpts).chhd(n_9)=V_hhd(:,n_9)'*V_hhd(:,n_9);
% Data(Cn_kpts).clhu(n_9)=V_lhu(:,n_9)'*V_lhu(:,n_9);
% Data(Cn_kpts).clhd(n_9)=V_lhd(:,n_9)'*V_lhd(:,n_9);
% Data(Cn_kpts).csou(n_9)=V_sou(:,n_9)'*V_sou(:,n_9);
% Data(Cn_kpts).csod(n_9)=V_sod(:,n_9)'*V_sod(:,n_9);
% %%%% Inner product |<Phi_c|Phi_v>|^2
% Data(Cn_kpts).Fc_Fv(n_8,n_9,1)=abs(V_ee(:,n_8)'*V_hhu(:,n_9))^2;
% Data(Cn_kpts).Fc_Fv(n_8,n_9,2)=abs(V_ee(:,n_8)'*V_hhd(:,n_9))^2;
% Data(Cn_kpts).Fc_Fv(n_8,n_9,3)=abs(V_ee(:,n_8)'*V_lhu(:,n_9))^2;
% Data(Cn_kpts).Fc_Fv(n_8,n_9,4)=abs(V_ee(:,n_8)'*V_lhd(:,n_9))^2;
% Data(Cn_kpts).Fc_Fv(n_8,n_9,5)=abs(V_ee(:,n_8)'*V_sou(:,n_9))^2;
% Data(Cn_kpts).Fc_Fv(n_8,n_9,6)=abs(V_ee(:,n_8)'*V_sod(:,n_9))^2;
% %%%% Calculate k=0 zone center eigenfunctions
% % if Cn_kpts == 1
% % V_HHu(:,n_9)=vectors(:,1);
% % V_HHd(:,n_9)=vectors(:,2);
% % V_LHu(:,n_9)=vectors(:,3);
% % V_LHd(:,n_9)=vectors(:,4);
% % V_SOu(:,n_9)=vectors(:,5);
% % V_SOd(:,n_9)=vectors(:,6);     
% % end
% end
% end
% end
Data(Cn_kpts).kx=kx;
Data(Cn_kpts).ky=ky;
eigs_time=eigs_time+toc;
waitbar((Cn_kpts/Tn_kpts),w);

% a=ful(Hmatrix_Holes);
end
end
end


close(w) 
disp('Total setup time'); disp(setup_time);
disp('Total eigs time'); disp(eigs_time);


%a=ful(Hmatrix_Holes);


%plot(

% figure(1)
% plot(1:21, (Eseperation(:,:,1)));
% figure(2)
% plot(1:21, (Eseperation(:,:,2)));
% figure(3)
% plot(1:21, (Eseperation(:,:,3)));
% figure(4)
% plot(1:21, (Eseperation(:,:,4)));
% figure(5)
% plot(1:21, (Eseperation(:,:,5)));
% figure(6)
% plot(1:21, (Eseperation(:,:,6)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check Valence band Dispersion
% Cn_kpts = 0;
% for y = 0:1:n_kpts-1
% for x = 0:1:n_kpts-1 
%     if x >= y
%     Cn_kpts = Cn_kpts +1;
%    
% %%%% Pull out CB and VB Dispersion energies in [01] and [11] directions
%     for w=1:n_Esub
%     if y == 0   
%     Ee_k_01(w,x+1)=Data(Cn_kpts).E_e(w);
%     Eh_k_01(w,x+1)=Data(Cn_kpts).E_h(w);
%     end
%     if x == y   
%     Ee_k_11(w,y+1)=Data(Cn_kpts).E_e(w);
%     Eh_k_11(w,y+1)=Data(Cn_kpts).E_h(w);
%     end
%     end
%     end
% end
% end
% k=0:(k_range/(n_kpts-1)):k_range;
% kpar01= ((2*pi)/a0)*sqrt(k.*k);
% kpar11= ((2*pi)/a0)*sqrt(k.*k + k.*k);
% 
% figure (1)
% hold on
% title('Valence Band Dispersion for 100A, GaAs-InGaAs Quantum Well')
% grid on
% plot(-kpar01, (Eh_k_01)','r', (kpar11), (Eh_k_11)','r');
% %axis([-0.06 0.06 -0.1 1.7]);
% text(-0.06,-0.01,'[01]');
% text(0.055,-0.01,'[11]');
% ylabel('Energy (eV)');
% xlabel('K|| (2 \pi /a_{0})'); %[A^0^-^1]
% hold off
% 
% figure (2)
% hold on
% title('Conduction Band Dispersion for 100A, SiGe-Ge Quantum Well')
% grid on
% plot(-kpar01, (Ee_k_01)','r', (kpar11), (Ee_k_11)','r');
% %axis([-0.06 0.06 -0.1 1.7]);
% text(-0.06,-0.01,'[01]');
% text(0.055,-0.01,'[11]');
% ylabel('Energy (eV)');
% xlabel('K|| (2 \pi /a_{0})'); %[A^0^-^1]
% hold off
% 
% figure (3)
% hold on
% title('Subband Dispersion for 100A, SiGe-Ge Quantum Well')
% grid on
% plot(-kpar01, (Eh_k_01)','b', (kpar11), (Eh_k_11)','b', -kpar01, (Ee_k_01)','r', (kpar11), (Ee_k_11)','r');
% %axis([-0.06 0.06 -0.1 1.7]);
% text(-0.06,-0.01,'[01]');
% text(0.055,-0.01,'[11]');
% ylabel('Energy (eV)');
% xlabel('K|| (2 \pi /a_{0})'); %[A^0^-^1]
% hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Stage 3.1 Save Data
% 
% if flag_1==0
% save(strcat(Sample_no,'.mat'),  ...
%         'a0',                   ... % Substrate lattice parameter
%         'n_Esub',               ... % Number of subbands
%         'k_range',              ... % Range of k_space covered
%         'n_kpts',               ... % Number of k points
%         'E_field',              ... % Applied electric field [V/um]
%         'E_p',                 ... % Optical matrix element
%         'n_r',                  ... % Refractive index of well region
%         'Data'                  ... % Data containing Eigen values and vectors
%      );
% 
% elseif flag_1 ==1
% save(strcat(Sample_no,'.mat'),  ...
%         'a0',                   ... % Substrate lattice parameter
%         'n_Esub',               ... % Number of subbands
%         'k_range',              ... % Range of k_space covered
%         'n_kpts',               ... % Number of k points
%         'E_field',              ... % Applied electric field [V/um]
%         'E_p',                 ... % Optical matrix element
%         'n_r',                  ... % Refractive index of well region
%         'Data',                 ... % Data containing Eigen values and vectors
%         'V_HHu',                ... %VB envelop fn, HH comp up
%         'V_HHd',                ... %VB envelop fn, HH comp down
%         'V_LHu',                ... %VB envelop fn, LH comp up
%         'V_LHd',                ... %VB envelop fn, LH comp down
%         'V_SOu',                ... %VB envelop fn, SO comp up
%         'V_SOd',                ... %VB envelop fn, SO comp down
%         'V_ee'                  ... %CB envelop fn
%                 );
% end