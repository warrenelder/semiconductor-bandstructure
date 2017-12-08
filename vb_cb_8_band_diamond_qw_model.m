%%%% Eigenstates and Eigenfunctions of electrons and holes in a single QW 
%%%% Author: Warren Elder

%%%% The double group exact envelope function equations are solved 
%%%% numerically, using a finite difference method and Brute force 
%%%% diagonalisation to obtain the eigenvalues and eigenfunctions of states
%%%% in a specified heterostructure
clear all 
clc
setup_time=0; % Counter for mesh set up
eigs_time=0;  % Counter for eigs diagonalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Control paramaters 
Sample_name = 'SiGe_QW';  % Sample name 
No_Bands    = 6;          % Number of bands; "6" 6-VB and 2-CB, "8" (CB,HH,LH,SO)  
substrate   = 'SiGe';     % Substrate material
Comp_vsub   = 0.90;       % Substrate composistion
n_layer     = 5;          % Number of epilayers
Material_layer1 = 'SiGe'; % Epilayer material; "SiGe" Si_{1-x}Ge_{x}
Material_layer2 = 'SiGe';
Material_layer3 = 'SiGe';
Material_layer4 = 'SiGe';
Material_layer5 = 'SiGe';
Material = [Material_layer1  % Define each eplilayer in heterostructure
            Material_layer2 
            Material_layer3 
            Material_layer4 
            Material_layer5];           
layer_comp  = [0.87 0.87 1 0.87 0.87]; % Composition of eplayers {x}
layer_npts  = [100 100 100 100 100];   % Number of mesh points in epilayer
layer_thck  = [100 100 100 100 100];   % Thickness (Angstrom) of epilayer
layer_cntr  = [0 0 1 0 0];  % Layer/Layers of interest (well) If there are multiple layers of interest, interest pointer must be seperated by a layer.
Strain_off  = 1;         % Strain perturbation; 0: off, 1: on
n_Esub      = 10;    	 % Number of eigenvalues or subbands to calculate
k_range    	= 0.06;   	 % Range of k_space covered as a percentage of 2pi/a
n_kpts      = 3;    	 % Number of k points in [hk] directions, excluding zone center
infboundary = 10;    	 % Choose boundary condition at end of heterostructure (eV)
Holes_Eval  = 0;         % Search for holes subband eigenvalues at energy (eV)
Elect_Eval  = 0.95;      % Search for electrons subband eigenvalues at energy (eV)
 
%%%% DO NOT MODIFY ANY CODE BEYOND THIS POINT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.1 Define scaling constants
% Define fundamental constants in S.I units. All energies are in eV from
% CODATA 2006 at http://physics.nist.gov/cuu/Constants/index.html
hbar   	= 1.054571628E-34;            	% Planck's constant [Js]
m0      = 9.10938215E-31;             	% Free electron mass[Kg]
e       = 1.602176487E-19;           	% Electron charge [C]
L       = 1E-10;                     	% Length Scale [Angstrom]
Const   = (hbar*hbar)/(2*m0*e*L*L);     % Scaling constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.2 Set up library and pre-allocate memory

L=repmat(struct('layer_num', 0,        ... %number of the layer
                'thick', 0,            ... %thickness of the layer
                'tot_thick', 0,        ... %total thickness of the region
                'layer_npts', 0,       ... %number of mesh points
                'tot_mesh_pts', 0,     ... %total mesh size
                'layer_mesh', 0,       ... %layer mesh sise
                'layer_cntr', 0,       ... %Layer of interest
                'comp', 0,             ... %layer composition
                'Me_g7m_cc', 0,        ... %Effective mass of electron
                'Me_g6m_cc', 0,        ... %Effective mass of electron
                'Z_g6m_vv', 0,         ... %Zeta paramter to G6- CB state
                'Z_g7m_vv', 0,         ... %Zeta paramter to G7- CB state
                'Z_g8m1_vv', 0,        ... %Zeta paramter to G8-1 CB state
                'Z_g8m2_vv', 0,        ... %Zeta paramter to G8-2 CB state
                'Z_g8m3_vv', 0,        ... %Zeta paramter to G8-3 CB state
                'Z_g7m_vs', 0,         ... %Zeta paramter to G7- CB state
                'Z_g8m1_vs', 0,        ... %Zeta paramter to G8-1 CB state
                'Z_g8m2_vs', 0,        ... %Zeta paramter to G8-2 CB state
                'Z_g7m_ss', 0,         ... %Zeta paramter to G7- CB state
                'Z_g8m_ss', 0,         ... %Zeta paramter to G8- CB state
                'Z_g6p_ccg6m', 0,      ... %Zeta paramter to G6- VB state
                'Z_g8p_ccg6m', 0,      ... %Zeta paramter to G8- VB state             
                'Z_g7p_ccg7m', 0,      ... %Zeta paramter to G7- VB state
                'Z_g8p_ccg7m', 0,      ... %Zeta paramter to G8- VB state  
                'Eg_unstrained', 0,    ... %Energy gap 
                'Delta', 0,            ... %spin orbit splitting
                'a',0,                 ... %Lattice parameter     
                'C11', 0,              ... %Elastic constant 'C11'
                'C12', 0,              ... %Elastic constant 'C12'                   
                'a_gap', 0,            ... %Deformation potential ('ac'-'av')
                'b', 0,                ... %'b' Deformation potential
                'offset_vb', 0,        ... %Valence band offset
                'offset_cb', 0,        ... %Conduction band offset
                'Luttbarp', 0,         ... %Luttbar
                'Mup', 0,              ... %Mu
                'deltap', 0,           ... %delta
                'Pip', 0,              ... %Pi
                'Sigmap', 0,           ... %Sigma
                'exx', 0,              ... %exx
                'eyy', 0,              ... %eyy
                'ezz', 0,              ... %ezz
                'Pe_hydro', 0,         ... %Pe_hydro
                'Qe_shear', 0,         ... %Qe_shear
                'h1', 0,               ... 
                'h2', 0,               ...
                'h3', 0,               ...
                'h4', 0,               ...
                'h5', 0,               ...
                'E_v1I', 0,            ...
                'E_v2I', 0,            ...
                'E_ssI', 0,            ...    
                'Pvv_cI', 0,           ...
                'Pvv_cI_con', 0,       ...
                'Qvv_cI', 0,           ...
                'Qvv_cI_con', 0,       ...
                'Svv_cI', 0,           ...
                'Svv_cI_dag', 0,       ...
                'Svv_cI_con', 0,       ...
                'Svv_cI_dag_con', 0,   ...
                'Rvv_cI', 0,           ...
                'Rvv_cI_dag', 0,       ...
                'Cvv_cI', 0,           ...
                'Cvv_cI_dag', 0,       ...
                'Zvv_cI', 0,           ...        
                'Zvv_cI_dag', 0,       ...
                'Qvs_cI', 0,           ...
                'Qvs_cI_dag', 0,       ...
                'Svs_cI', 0,           ...
                'Svs_cI_dag', 0,       ...
                'Svs_cI_con', 0,       ...
                'Svs_cI_dag_con', 0,   ...
                'SIGvs_cI', 0,         ...
                'SIGvs_cI_dag', 0,     ...
                'SIGvs_cI_con', 0,     ...
                'SIGvs_cI_dag_con', 0, ...
                'Rvs_cI', 0,           ...
                'Rvs_cI_dag', 0,       ...
                'Rvs_cI_con', 0,       ...
                'Rvs_cI_dag_con', 0,   ...
                'Ps_cI', 0,            ...
                'Ps_cI_con', 0,        ...
                'Cs_cI', 0,            ...
                'Cs_cI_dag', 0,        ...
                'Pvv_pI', 0,           ...
                'Pvv_pI_con', 0,       ...
                'Qvv_pI', 0,           ...
                'Qvv_pI_con', 0,       ...
                'Svv_pI', 0,           ...
                'Svv_pI_dag', 0,       ...
                'Svv_pI_con', 0,       ...
                'Svv_pI_dag_con', 0,   ...
                'Rvv_pI', 0,           ...
                'Rvv_pI_dag', 0,       ...
                'Cvv_pI', 0,           ...
                'Cvv_pI_dag', 0,       ...
                'Zvv_pI', 0,           ...        
                'Zvv_pI_dag', 0,       ...
                'Qvs_pI', 0,           ...
                'Qvs_pI_dag_con', 0,   ...
                'Qvs_pI_con', 0,       ...
                'Qvs_pI_dag', 0,       ...
                'Svs_pI', 0,           ...
                'Svs_pI_dag', 0,       ...
                'Svs_pI_con', 0,       ...
                'Svs_pI_dag_con', 0,   ...
                'SIGvs_pI', 0,         ...
                'SIGvs_pI_dag', 0,     ...
                'SIGvs_pI_con', 0,     ...
                'SIGvs_pI_dag_con', 0, ...
                'Rvs_pI', 0,           ...
                'Rvs_pI_dag', 0,       ...
                'Rvs_pI_con', 0,       ...
                'Rvs_pI_dag_con', 0,   ...
                'Ps_pI', 0,            ...
                'Ps_pI_con', 0,        ...
                'Cs_pI', 0,            ...
                'Cs_pI_dag', 0,        ...
                'Pvv_mI', 0,           ...
                'Pvv_mI_con', 0,       ...
                'Qvv_mI', 0,           ...
                'Qvv_mI_con', 0,       ...
                'Svv_mI', 0,           ...
                'Svv_mI_dag', 0,       ...
                'Svv_mI_con', 0,       ...
                'Svv_mI_dag_con', 0,   ...
                'Rvv_mI', 0,           ...
                'Rvv_mI_dag', 0,       ...
                'Cvv_mI', 0,           ...
                'Cvv_mI_dag', 0,       ...
                'Zvv_mI', 0,           ...        
                'Zvv_mI_dag', 0,       ...
                'Qvs_mI', 0,           ...
                'Qvs_mI_dag_con', 0,   ...
                'Qvs_mI_con', 0,       ...
                'Qvs_mI_dag', 0,       ...
                'Svs_mI', 0,           ...
                'Svs_mI_dag', 0,       ...
                'Svs_mI_con', 0,       ...
                'Svs_mI_dag_con', 0,   ...
                'SIGvs_mI', 0,         ...
                'SIGvs_mI_dag', 0,     ...
                'SIGvs_mI_con', 0,     ...
                'SIGvs_mI_dag_con', 0, ...
                'Rvs_mI', 0,           ...
                'Rvs_mI_dag', 0,       ...
                'Rvs_mI_con', 0,       ...
                'Rvs_mI_dag_con', 0,   ...
                'Ps_mI', 0,            ...
                'Ps_mI_con', 0,        ...
                'Cs_mI', 0,            ...
                'Cs_mI_dag', 0,        ...
                'E_ccIL', 0,           ...
                'E_v1IL', 0,           ...
                'E_v2IL', 0,           ...
                'E_ssIL', 0,           ...
                'Pvv_cIL', 0,          ...
                'Pvv_cIL_con', 0,      ...
                'Qvv_cIL', 0,          ...
                'Qvv_cIL_con', 0,      ...
                'Svv_cIL', 0,          ...
                'Svv_cIL_dag', 0,      ...
                'Svv_cIL_con', 0,      ...
                'Svv_cIL_dag_con', 0,  ...
                'Rvv_cIL', 0,          ...
                'Rvv_cIL_dag', 0,      ...
                'Cvv_cIL', 0,          ...
                'Cvv_cIL_dag', 0,      ...
                'Zvv_cIL', 0,          ...    
                'Zvv_cIL_dag', 0,      ...
                'Qvs_cIL', 0,          ...
                'Qvs_cIL_dag_con', 0,  ...
                'Qvs_cIL_con', 0,      ...
                'Qvs_cIL_dag', 0,      ...
                'Svs_cIL', 0,          ...
                'Svs_cIL_dag', 0,      ...
                'Svs_cIL_con', 0,      ...
                'Svs_cIL_dag_con', 0,  ...
                'SIGvs_cIL', 0,        ...
                'SIGvs_cIL_dag', 0,    ...
                'SIGvs_cIL_con', 0,    ...
                'SIGvs_cIL_dag_con', 0, ...
                'Rvs_cIL', 0,          ...
                'Rvs_cIL_dag', 0,      ...
                'Rvs_cIL_con', 0,      ...
                'Rvs_cIL_dag_con', 0,  ...
                'Ps_cIL', 0,           ...
                'Ps_cIL_con', 0,       ...
                'Cs_cIL', 0,           ...
                'Cs_cIL_dag', 0,       ...
                'E_ccIR', 0,           ...
                'E_v1IR', 0,           ...
                'E_v2IR', 0,           ...
                'E_ssIR', 0,           ...
                'Pvv_cIR', 0,          ...
                'Pvv_cIR_con', 0,      ...
                'Qvv_cIR', 0,          ...
                'Qvv_cIR_con', 0,      ...
                'Svv_cIR', 0,          ...
                'Svv_cIR_dag', 0,      ...
                'Svv_cIR_con', 0,      ...
                'Svv_cIR_dag_con', 0,  ...
                'Rvv_cIR', 0,          ...
                'Rvv_cIR_dag', 0,      ...
                'Cvv_cIR', 0,          ...
                'Cvv_cIR_dag', 0,      ...
                'Zvv_cIR', 0,          ...
                'Zvv_cIR_dag', 0,      ...
                'Qvs_cIR', 0,          ...
                'Qvs_cIR_dag_con', 0,  ...
                'Qvs_cIR_con', 0,      ...
                'Qvs_cIR_dag', 0,      ...
                'Svs_cIR', 0,          ...
                'Svs_cIR_dag', 0,      ...
                'Svs_cIR_con', 0,      ...
                'Svs_cIR_dag_con', 0,  ...
                'SIGvs_cIR', 0,        ...
                'SIGvs_cIR_dag', 0,    ...
                'SIGvs_cIR_con', 0,    ...
                'SIGvs_cIR_dag_con', 0, ...
                'Rvs_cIR', 0,          ...
                'Rvs_cIR_dag', 0,      ...
                'Rvs_cIR_con', 0,      ...
                'Rvs_cIR_dag_con', 0,  ...
                'Ps_cIR', 0,           ...
                'Ps_cIR_con', 0,       ...
                'Cs_cIR', 0,           ...
                'Cs_cIR_dag', 0,       ...  
                'E_v1J', 0,            ...
                'E_v2J', 0,            ...
                'E_ssJ', 0,            ...    
                'Pvv_cJ', 0,           ...
                'Pvv_cJ_con', 0,       ...
                'Qvv_cJ', 0,           ...
                'Qvv_cJ_con', 0,       ...
                'Svv_cJ', 0,           ...
                'Svv_cJ_dag', 0,       ...
                'Svv_cJ_con', 0,       ...
                'Svv_cJ_dag_con', 0,   ...
                'Rvv_cJ', 0,           ...
                'Rvv_cJ_dag', 0,       ...
                'Cvv_cJ', 0,           ...
                'Cvv_cJ_dag', 0,       ...
                'Zvv_cJ', 0,           ...        
                'Zvv_cJ_dag', 0,       ...
                'Qvs_cJ', 0,           ...
                'Qvs_cJ_dag', 0,       ...
                'Svs_cJ', 0,           ...
                'Svs_cJ_dag', 0,       ...
                'Svs_cJ_con', 0,       ...
                'Svs_cJ_dag_con', 0,   ...
                'SIGvs_cJ', 0,         ...
                'SIGvs_cJ_dag', 0,     ...
                'SIGvs_cJ_con', 0,     ...
                'SIGvs_cJ_dag_con', 0, ...
                'Rvs_cJ', 0,           ...
                'Rvs_cJ_dag', 0,       ...
                'Rvs_cJ_con', 0,       ...
                'Rvs_cJ_dag_con', 0,   ...
                'Ps_cJ', 0,            ...
                'Ps_cJ_con', 0,        ...
                'Cs_cJ', 0,            ...
                'Cs_cJ_dag', 0,        ...
                'Pvv_pJ', 0,           ...
                'Pvv_pJ_con', 0,       ...
                'Qvv_pJ', 0,           ...
                'Qvv_pJ_con', 0,       ...
                'Svv_pJ', 0,           ...
                'Svv_pJ_dag', 0,       ...
                'Svv_pJ_con', 0,       ...
                'Svv_pJ_dag_con', 0,   ...
                'Rvv_pJ', 0,           ...
                'Rvv_pJ_dag', 0,       ...
                'Cvv_pJ', 0,           ...
                'Cvv_pJ_dag', 0,       ...
                'Zvv_pJ', 0,           ...        
                'Zvv_pJ_dag', 0,       ...
                'Qvs_pJ', 0,           ...
                'Qvs_pJ_dag_con', 0,   ...
                'Qvs_pJ_con', 0,       ...
                'Qvs_pJ_dag', 0,       ...
                'Svs_pJ', 0,           ...
                'Svs_pJ_dag', 0,       ...
                'Svs_pJ_con', 0,       ...
                'Svs_pJ_dag_con', 0,   ...
                'SIGvs_pJ', 0,         ...
                'SIGvs_pJ_dag', 0,     ...
                'SIGvs_pJ_con', 0,     ...
                'SIGvs_pJ_dag_con', 0, ...
                'Rvs_pJ', 0,           ...
                'Rvs_pJ_dag', 0,       ...
                'Rvs_pJ_con', 0,       ...
                'Rvs_pJ_dag_con', 0,   ...
                'Ps_pJ', 0,            ...
                'Ps_pJ_con', 0,        ...
                'Cs_pJ', 0,            ...
                'Cs_pJ_dag', 0,        ...
                'Pvv_mJ', 0,           ...
                'Pvv_mJ_con', 0,       ...
                'Qvv_mJ', 0,           ...
                'Qvv_mJ_con', 0,       ...
                'Svv_mJ', 0,           ...
                'Svv_mJ_dag', 0,       ...
                'Svv_mj_con', 0,       ...
                'Svv_mJ_dag_con', 0,   ...
                'Rvv_mJ', 0,           ...
                'Rvv_mJ_dag', 0,       ...
                'Cvv_mJ', 0,           ...
                'Cvv_mJ_dag', 0,       ...
                'Zvv_mJ', 0,           ...        
                'Zvv_mJ_dag', 0,       ...
                'Qvs_mJ', 0,           ...
                'Qvs_mJ_dag_con', 0,   ...
                'Qvs_mJ_con', 0,       ...
                'Qvs_mJ_dag', 0,       ...
                'Svs_mJ', 0,           ...
                'Svs_mJ_dag', 0,       ...
                'Svs_mJ_con', 0,       ...
                'Svs_mJ_dag_con', 0,   ...
                'SIGvs_mJ', 0,         ...
                'SIGvs_mJ_dag', 0,     ...
                'SIGvs_mJ_con', 0,     ...
                'SIGvs_mJ_dag_con', 0, ...
                'Rvs_mJ', 0,           ...
                'Rvs_mJ_dag', 0,       ...
                'Rvs_mJ_con', 0,       ...
                'Rvs_mJ_dag_con', 0,   ...
                'Ps_mJ', 0,            ...
                'Ps_mJ_con', 0,        ...
                'Cs_mJ', 0,            ...
                'Cs_mJ_dag', 0         ... 
                                          ), 1, n_layer);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.3 Define material parameters and number of bands in model
layer_num = 1 : n_layer;
for n_1 = 1 : n_layer
L(n_1).layer_num  = layer_num(n_1); 
L(n_1).layer_npts = layer_npts(n_1);          
L(n_1).comp       = layer_comp(n_1);
tot_mesh_pts      = sum(layer_npts); 
tot_thick         = sum(layer_thck); 
L(n_1).layer_mesh = layer_thck(n_1)./layer_npts(n_1);
L(n_1).layer_cntr = layer_cntr(n_1);
L(n_1).layer_thck = layer_thck(n_1);
L(n_1).h          = L(n_1).layer_mesh;

switch Material(n_1,:)
%%%% Material parameters for different material layers
case 'SiGe'

a0 = 5.431 + 0.1992*Comp_vsub +0.02733*Comp_vsub*Comp_vsub;

L(n_1).Xi_G8p_G6m   =  (1.118313462317252 + 0.082885938281999*layer_comp(n_1));
L(n_1).Xi_G8p_G7m   =  (2.505261703985227 - 0.253783087969172*layer_comp(n_1));
L(n_1).Xi_G8p_G8m_1 =  (3.177262973063451 + 0.394834450020080*layer_comp(n_1));
L(n_1).Xi_G8p_G8m_2 = -(0.727680986725345 + 0.172000274182098*layer_comp(n_1));
L(n_1).Xi_G7p_G7m   = -(0.660000000000000 + 1.014000000000000*layer_comp(n_1));
L(n_1).Xi_G7p_G8m   = -(1.937000000000000 + 0.888000000000000*layer_comp(n_1));

L(n_1).E_G8p = (0);
L(n_1).E_G7p = (-0.044 - 0.253*layer_comp(n_1));
L(n_1).E_G7m = (4.185 - 3.298*layer_comp(n_1));
L(n_1).E_G6m = (3.335 - 0.329*layer_comp(n_1));
L(n_1).E_G8m = (3.365 - 0.159*layer_comp(n_1));

if layer_comp(n_1)<=0.12 
L(n_1).Eg_unstrained = L(n_1).E_G6m -L(n_1).E_G8p;
elseif layer_comp(n_1)>0.12
L(n_1).Eg_unstrained = L(n_1).E_G7m -L(n_1).E_G8p;    
end
L(n_1).Delta = -(-0.044 - 0.253*layer_comp(n_1));  
L(n_1).a     = (5.431 + 0.1992*layer_comp(n_1)+0.02733*layer_comp(n_1)*layer_comp(n_1));
L(n_1).C11   = (16.577 - 3.724*layer_comp(n_1));
L(n_1).C12   = (6.393 - 1.565*layer_comp(n_1));
L(n_1).C44   = (7.962 - 6.680*layer_comp(n_1));
L(n_1).ac    = (-10.06 + 2.23*layer_comp(n_1));
L(n_1).av    = (2.38 - 0.15*layer_comp(n_1)); 
L(n_1).a_gap = (-12.19 + 0.54*layer_comp(n_1));
L(n_1).bv    = -(-2.10 - 0.76*layer_comp(n_1));
L(n_1).dv    = (-4.85 - 0.43*layer_comp(n_1));

end 
%%%% Second order interaction paramters and Luttinger paramters for 6- and
%%%% 8-band models
if No_Bands == 6
L(n_1).Z_g6m_vv  = (((L(n_1).Xi_G8p_G6m)^2)/(L(n_1).E_G6m)); 
L(n_1).Z_g7m_vv  = (((L(n_1).Xi_G8p_G7m)^2)/(L(n_1).E_G7m)); 
L(n_1).Z_g8m1_vv = (((L(n_1).Xi_G8p_G8m_1)^2)/(L(n_1).E_G8m));
L(n_1).Z_g8m2_vv = (((L(n_1).Xi_G8p_G8m_2)^2)/(L(n_1).E_G8m));
L(n_1).Z_g8m3_vv = (((L(n_1).Xi_G8p_G8m_1)*(L(n_1).Xi_G8p_G8m_2))/(L(n_1).E_G8m));
L(n_1).Z_g7m_vs  = (((1/2)*(L(n_1).Xi_G8p_G7m)*(L(n_1).Xi_G7p_G7m))*((1/(L(n_1).E_G7m)) + (1/(L(n_1).E_G7m-L(n_1).E_G7p))));
L(n_1).Z_g8m1_vs = (((1/2)*(L(n_1).Xi_G8p_G8m_1)*(L(n_1).Xi_G7p_G8m))*((1/(L(n_1).E_G8m)) + (1/(L(n_1).E_G8m-L(n_1).E_G7p))));
L(n_1).Z_g8m2_vs = (((1/2)*(L(n_1).Xi_G8p_G8m_2)*(L(n_1).Xi_G7p_G8m))*((1/(L(n_1).E_G8m)) + (1/(L(n_1).E_G8m-L(n_1).E_G7p))));
L(n_1).Z_g7m_ss  = (((L(n_1).Xi_G7p_G7m)^2)/((L(n_1).E_G7m-L(n_1).E_G7p)));
L(n_1).Z_g8m_ss  = (((L(n_1).Xi_G7p_G8m)^2)/((L(n_1).E_G8m-L(n_1).E_G7p)));
L(n_1).Z_g7p_cc  = ((L(n_1).Xi_G7p_G7m)^2)/(-(L(n_1).E_G7m-L(n_1).E_G7p));
L(n_1).Z_g8p_cc  = ((L(n_1).Xi_G8p_G7m)^2)/(-L(n_1).E_G7m);
L(n_1).Xi_cv     = 0;
L(n_1).Xi_cs     = 0;
elseif No_Bands == 8
L(n_1).Z_g6m_vv  = (((L(n_1).Xi_G8p_G6m)^2)/(L(n_1).E_G6m)); 
L(n_1).Z_g7m_vv  = 0; 
L(n_1).Z_g8m1_vv = (((L(n_1).Xi_G8p_G8m_1)^2)/(L(n_1).E_G8m));
L(n_1).Z_g8m2_vv = (((L(n_1).Xi_G8p_G8m_2)^2)/(L(n_1).E_G8m));
L(n_1).Z_g8m3_vv = (((L(n_1).Xi_G8p_G8m_1)*(L(n_1).Xi_G8p_G8m_2))/(L(n_1).E_G8m));
L(n_1).Z_g7m_vs  = 0;
L(n_1).Z_g8m1_vs = (((1/2)*(L(n_1).Xi_G8p_G8m_1)*(L(n_1).Xi_G7p_G8m))*((1/(L(n_1).E_G8m)) + (1/(L(n_1).E_G8m-L(n_1).E_G7p))));
L(n_1).Z_g8m2_vs = (((1/2)*(L(n_1).Xi_G8p_G8m_2)*(L(n_1).Xi_G7p_G8m))*((1/(L(n_1).E_G8m)) + (1/(L(n_1).E_G8m-L(n_1).E_G7p))));
L(n_1).Z_g7m_ss  = 0;
L(n_1).Z_g8m_ss  = (((L(n_1).Xi_G7p_G8m)^2)/((L(n_1).E_G8m-L(n_1).E_G7p)));
L(n_1).Z_g7p_cc  = 0;
L(n_1).Z_g8p_cc  = 0;
L(n_1).Xi_cv     = (L(n_1).Xi_G8p_G7m);
L(n_1).Xi_cs     = (L(n_1).Xi_G7p_G7m);    
end
L(n_1).L1vv      = (2*L(n_1).Z_g6m_vv + 2*L(n_1).Z_g7m_vv  +  L(n_1).Z_g8m1_vv  + 8*L(n_1).Z_g8m2_vv + 4*L(n_1).Z_g8m3_vv-1); 
L(n_1).L2vv      = (L(n_1).Z_g6m_vv - L(n_1).Z_g7m_vv  -  4*L(n_1).Z_g8m2_vv - 2*L(n_1).Z_g8m3_vv);
L(n_1).L3vv      = (L(n_1).Z_g6m_vv + L(n_1).Z_g7m_vv  -  2*L(n_1).Z_g8m2_vv);
L(n_1).L2vs      = (sqrt(1/2)*(L(n_1).Z_g7m_vs - L(n_1).Z_g8m1_vs - 4*L(n_1).Z_g8m2_vs));
L(n_1).L3vs      = (sqrt(1/2)*(L(n_1).Z_g7m_vs + L(n_1).Z_g8m1_vs + 2*L(n_1).Z_g8m2_vs)); 
L(n_1).L1ss      = (L(n_1).Z_g7m_ss + 4*L(n_1).Z_g8m_ss -1); 
L(n_1).L1eg7m    = (L(n_1).Z_g7p_cc + 4*L(n_1).Z_g8p_cc -1);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 1.4 Call material parameters, set material mesh, set data file
%%%% and pre-allocate memory
n_vector        = tot_mesh_pts*8;
order_ptr       = zeros(n_Esub,1);
Tn_kpts         = n_kpts*(n_kpts+1)/2;
Cn_kpts         = 0;
I1              = diag([ones(2,1) ; -ones(6,1)]);
I2              = diag([zeros(2,1) ; -zeros(6,1)]);
infboundary_h1  = infboundary*I1;
infboundary_h2  = infboundary*I2;
E_valelec1      = ones(n_Esub,n_kpts);
E_valholes1     = ones(n_Esub,n_kpts);
E_valCB{n_kpts} = zeros(n_Esub,n_kpts);
E_valVB{n_kpts} = zeros(n_Esub,n_kpts);
E_H             = zeros(n_kpts,n_kpts);   
E_E             = zeros(n_kpts,n_kpts);
V_ecu_h         = zeros((tot_mesh_pts+1),n_Esub);
V_ecd_h         = zeros((tot_mesh_pts+1),n_Esub);
V_hhu_h         = zeros((tot_mesh_pts+1),n_Esub);
V_hhd_h         = zeros((tot_mesh_pts+1),n_Esub);
V_lhu_h         = zeros((tot_mesh_pts+1),n_Esub);
V_lhd_h         = zeros((tot_mesh_pts+1),n_Esub);
V_sou_h         = zeros((tot_mesh_pts+1),n_Esub);
V_sod_h         = zeros((tot_mesh_pts+1),n_Esub);
V_ecu_e         = zeros((tot_mesh_pts+1),n_Esub);
V_ecd_e         = zeros((tot_mesh_pts+1),n_Esub);
V_hhu_e         = zeros((tot_mesh_pts+1),n_Esub);
V_hhd_e         = zeros((tot_mesh_pts+1),n_Esub);
V_lhu_e         = zeros((tot_mesh_pts+1),n_Esub);
V_lhd_e         = zeros((tot_mesh_pts+1),n_Esub);
V_sou_e         = zeros((tot_mesh_pts+1),n_Esub);
V_sod_e         = zeros((tot_mesh_pts+1),n_Esub);
%%%% Defines the calculated Data Structure

Data=repmat(struct( 'kx', 0,                   ... %x comp of wave vector
                    'ky', 0,                   ... %y comp of wave vector
                    'E_H', zeros(1,n_Esub),    ... %eigen energy electrons
                    'E_E', zeros(1,n_Esub),    ... %eigen energy holes
                    'V_hhu', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, HH comp up
                    'V_hhd', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, HH comp down
                    'V_lhu', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, LH comp up
                    'V_lhd', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, LH comp down
                    'V_sou', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, SO comp up
                    'V_sod', zeros(tot_mesh_pts+1, n_Esub),    ... %VB envelop fn, SO comp down
                    'V_ee', zeros(tot_mesh_pts+1, n_Esub)     ... %CB envelop fn
                            ), 1, Tn_kpts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.1 Calculate the Eigenenergies and Eigenfunctions at different 
%%%% K points from X to L crystallagraphic directions
w = waitbar(0,'Wait wait wait wait...'); % Implement waitbar
for y = 0:1:n_kpts-1
for x = 0:1:n_kpts-1     
if (x>=y)
tic;
%%%% Stage 2.2 Define the in-plane wavevector components and scale with the
%%%% lattice constant of the substrate
Cn_kpts = Cn_kpts+1;
k_x  = x*(k_range/(n_kpts-1));
k_y  = y*(k_range/(n_kpts-1));
kx   = (2*pi/a0)*k_x;
ky   = (2*pi/a0)*k_y;
kpar = (2*pi/a0)*sqrt(k_x*k_x + k_y*k_y);
k_p  = (2*pi/a0)*(k_x + 1i*k_y);
k_m  = (2*pi/a0)*(k_x - 1i*k_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.2  Set model solid theory calculated band edge and define 
%%%% strain conditions 
for n_2 =1:n_layer
L(n_2).epar  = (a0-L(n_2).a)/L(n_2).a;
L(n_2).eperp = ((-2*L(n_2).C12)/L(n_2).C11)*((a0-L(n_2).a)/L(n_2).a);
Ev_Si        = 0;
L(n_2).Ev_av = Ev_Si+((0.47 + 0.06*Comp_vsub)*(layer_comp(n_2)))...
+L(n_2).av*(3*L(n_2).epar-(2*L(n_2).epar +L(n_2).eperp))+(1/3)*L(n_2).Delta;
L(n_2).E_v1  = L(n_2).Ev_av + L(n_2).E_G8p;
L(n_2).E_v2  = L(n_2).Ev_av + L(n_2).E_G8p;
L(n_2).E_ss  = L(n_2).Ev_av + L(n_2).E_G7p;
L(n_2).E_cc  = L(n_2).Ev_av + L(n_2).E_G7m + (L(n_2).ac-L(n_2).av)*(3*L(n_2).epar);      
if Strain_off == 0
L(n_2).euu = 0;
L(n_2).evv = 0;
L(n_2).eww = 0;
L(n_2).euv = 0;
L(n_2).evw = 0;
L(n_2).ewu = 0; 
else
L(n_2).euu = (L(n_2).epar);
L(n_2).evv = (L(n_2).epar);
L(n_2).eww = (L(n_2).eperp);
L(n_2).euv = 0;
L(n_2).evw = 0;
L(n_2).ewu = 0;
end
L(n_2).E_p = (L(n_2).av*(L(n_2).euu + L(n_2).evv + L(n_2).eww - 3*L(n_2).epar));
L(n_2).E_q = ((1/2)*L(n_2).bv*(L(n_2).euu + L(n_2).evv - 2*L(n_2).eww));
L(n_2).E_s = (L(n_2).dv*(L(n_2).ewu - 1i*L(n_2).evw));
L(n_2).E_r = (((-sqrt(3)/2)*L(n_2).bv*(L(n_2).euu - L(n_2).evv) + 1i*L(n_2).dv*L(n_2).euv));
L(n_2).E_c = (L(n_2).ac*(L(n_2).euu + L(n_2).evv + L(n_2).eww - 3*L(n_2).epar));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.3 Set Conduction and Valence band profiles at k=0
for n_3 = 1 : n_layer        
L(n_3).B_edge =  [(L(n_3).E_cc)
                  (L(n_3).E_cc)
                  (L(n_3).E_v1) 
                  (L(n_3).E_v2) 
                  (L(n_3).E_v2)  
                  (L(n_3).E_v1)   
                  (L(n_3).E_ss)     
                  (L(n_3).E_ss)];  
L(n_3).V_pot =  [(L(n_3).E_cc + L(n_3).E_c) (0)                        (0)                                     (0)                                        (0)                                      (0)                                        (0)                           (0)                        
                 (0)                        (L(n_3).E_cc + L(n_3).E_c) (0)                                     (0)                                        (0)                                      (0)                                        (0)                           (0) 
                 (0)                        (0)                        (L(n_3).E_v1 + L(n_3).E_p + L(n_3).E_q) (-L(n_3).E_s)                              (L(n_3).E_r)                             (0)                                        (sqrt(3/2)*L(n_3).E_s')       (sqrt(2)*L(n_3).E_q) 
                 (0)                        (0)                        (-L(n_3).E_s')                          (L(n_3).E_v2 + L(n_3).E_p - L(n_3).E_q) 	  (0)                                      (L(n_3).E_r)                               (sqrt(2)*L(n_3).E_r')         (sqrt(1/2)*L(n_3).E_s)
                 (0)                        (0)                        (L(n_3).E_r')                           (0)                                        (L(n_3).E_v2 + L(n_3).E_p - L(n_3).E_q)  (L(n_3).E_s)                               (sqrt(1/2)*L(n_3).E_s)        (-sqrt(2)*L(n_3).E_r)  
                 (0)                        (0)                        (0)                                     (L(n_3).E_r')                              (L(n_3).E_s')                            (L(n_3).E_v1 + L(n_3).E_p + L(n_3).E_q)    (-sqrt(2)*L(n_3).E_q)         (sqrt(3/2)*L(n_3).E_s)
                 (0)                        (0)                        (sqrt(3/2)*L(n_3).E_s)                  (sqrt(2)*L(n_3).E_r)                       (sqrt(1/2)*L(n_3).E_s')                  (-sqrt(2)*L(n_3).E_q)                      (L(n_3).E_ss + L(n_3).E_p)    (0)     
                 (0)                        (0)                        (sqrt(2)*L(n_3).E_q')                   (sqrt(1/2)*L(n_3).E_s)                     (-sqrt(2)*L(n_3).E_r')                   (sqrt(3/2)*L(n_3).E_s')                    (0)                           (L(n_3).E_ss + L(n_3).E_p)];
end 
if Cn_kpts == 1
for u_3 = 1:n_layer
u_4 = layer_cntr(1,u_3);
if u_4>0
Layermeshpts=L(u_3).layer_npts;

B_P(u_3).Band_profile_0 =  L(u_3).B_edge;
B_P(u_3).Band_profile_1 =  (eig(L(u_3).V_pot));  

B_P(u_3).E_celectron = B_P(u_3).Band_profile_1(7)*ones(1,Layermeshpts); 
B_P(u_3).E_heavyhole = B_P(u_3).Band_profile_1(5)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole = B_P(u_3).Band_profile_1(3)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit = B_P(u_3).Band_profile_1(1)*ones(1,Layermeshpts);

B_P(u_3).E_celectron_unp = B_P(u_3).Band_profile_0(1)*ones(1,Layermeshpts); 
B_P(u_3).E_heavyhole_unp = B_P(u_3).Band_profile_0(3)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole_unp = B_P(u_3).Band_profile_0(5)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit_unp = B_P(u_3).Band_profile_0(7)*ones(1,Layermeshpts);
else
Layermeshpts=L(u_3).layer_npts;

B_P(u_3).Band_profile_0 =  L(u_3).B_edge;
B_P(u_3).Band_profile_1 =  (eig(L(u_3).V_pot)); 

B_P(u_3).E_celectron = B_P(u_3).Band_profile_1(7)*ones(1,Layermeshpts);
B_P(u_3).E_heavyhole = B_P(u_3).Band_profile_1(3)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole = B_P(u_3).Band_profile_1(5)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit = B_P(u_3).Band_profile_1(1)*ones(1,Layermeshpts);

B_P(u_3).E_celectron_unp = B_P(u_3).Band_profile_0(1)*ones(1,Layermeshpts);
B_P(u_3).E_heavyhole_unp = B_P(u_3).Band_profile_0(3)*ones(1,Layermeshpts);                       
B_P(u_3).E_lighthole_unp = B_P(u_3).Band_profile_0(5)*ones(1,Layermeshpts);
B_P(u_3).E_spinorbit_unp = B_P(u_3).Band_profile_0(7)*ones(1,Layermeshpts);
end
if L(u_3).layer_cntr == 1   
Offset =  B_P(u_3).E_heavyhole(1,1);   
end
end
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.4 This loop calculate the Block Tridiagonal Matrix Hamiltonian 
%%%% elements. The forward-backward (center difference finite difference 
%%%% approximation is used to obtain discretitized values of the 
%%%% Hamiltonian at each of the N mesh points. These values are then 
%%%% re-cast into a cell array forming a block tridiagonal matrix 
%%%% consisting of cells for diagonal-, sub- and super-diagonal blocks at 
%%%% each mesh point. Below provides examples of treatment of different 
%%%% operators under finite difference scheme, where "gi" is the material 
%%%% dependent parameter and his the distance between mesh points.

% L(n_4).Aij_cI  = (k^2)*(((gi) + (gi-1))/2) -(-2)*(((gi)*L(n_4-1).h + (gi-1)*L(n_4).h)/L(n_4).h1); 
% L(n_4).Aij_cIR = (k^2)*((gi)/2)            -(-2)*(((3/2)*(gi)   + (1/4)*((gi) + (gi-1)))/L(n_4).h3); 
% L(n_4).Aij_cIL = (k^2)*((gi-1)/2)          -(-2)*(((3/2)*(gi-1) + (1/4)*((gi) + (gi-1)))/L(n_4).h3); 
% L(n_4).Aij_pI  =                           +(-2)*((((gi)   + ((gi-1) - (gi))/4)*L(n_4).h)/L(n_4).h1)      + ((gi) + (gi-1))/(2*L(n_4).h2) + ((gi)/(L(n_4).h2)); 
% L(n_4).Aij_mI  =                           +(-2)*((((gi-1) - ((gi-1) -  (gi))/4)*L(n_4-1).h)/L(n_4).h1)    - ((gi) + (gi-1))/(2*L(n_4).h2) - ((gi-1)/(L(n_4).h2)); 
for n_4 = 1 : n_layer
if (n_4 == 1) 
%%%% Set up first interface matrix at start of heterostructure
%%%% Stage 1 Center Fj=i points 
L(n_4).E_ccI = infboundary + Offset;
L(n_4).E_v1I = -infboundary + Offset;
L(n_4).E_v2I = -infboundary + Offset;
L(n_4).E_ssI = -infboundary + Offset;  
L(n_4).E_cI         = 0;
L(n_4).E_pI         = 0;
L(n_4).E_qI         = 0;
L(n_4).E_sI         = 0;
L(n_4).E_rI         = 0;
L(n_4).Pe_cI        = 0;
L(n_4).Pe_cI_con    = 0;
L(n_4).Ce_cI        = 0;
L(n_4).Ce_cI_dag    = 0;
L(n_4).Te_cI       = 0;
L(n_4).Te_cI_dag   = 0;
L(n_4).Ue_cI       = 0;
L(n_4).We_cI       = 0;
L(n_4).Ve_cI       = 0;
L(n_4).Ve_cI_dag   = 0;   
L(n_4).Pvv_cI       = 0;
L(n_4).Pvv_cI_con   = 0;
L(n_4).Qvv_cI       = 0;
L(n_4).Qvv_cI_con   = 0;
L(n_4).Svv_cI       = 0;
L(n_4).Svv_cI_dag   = 0;
L(n_4).Svv_cI_con   = 0;
L(n_4).Svv_cI_dag_con = 0;
L(n_4).Rvv_cI       = 0;
L(n_4).Rvv_cI_dag   = 0;
L(n_4).Cvv_cI       = 0;
L(n_4).Cvv_cI_dag   = 0;
L(n_4).Zvv_cI       = 0;        
L(n_4).Zvv_cI_dag   = 0;
L(n_4).Qvs_cI       = 0;
L(n_4).Qvs_cI_dag   = 0;
L(n_4).Svs_cI       = 0;
L(n_4).Svs_cI_dag   = 0;
L(n_4).Svs_cI_con   = 0;
L(n_4).Svs_cI_dag_con = 0;
L(n_4).SIGvs_cI     = 0;
L(n_4).SIGvs_cI_dag = 0;
L(n_4).SIGvs_cI_con = 0;
L(n_4).SIGvs_cI_dag_con = 0;
L(n_4).Rvs_cI       = 0;
L(n_4).Rvs_cI_dag   = 0;
L(n_4).Rvs_cI_con   = 0;
L(n_4).Rvs_cI_dag_con = 0;
L(n_4).Ps_cI        = 0;
L(n_4).Ps_cI_con    = 0;
L(n_4).Cs_cI        = 0;
L(n_4).Cs_cI_dag    = 0;
%%%% Stage 2 Forward F{j=i+1) points
L(n_4).Pe_pI        = 0;
L(n_4).Pe_pI_con    = 0;
L(n_4).Ce_pI        = 0;
L(n_4).Ce_pI_dag    = 0;
L(n_4).Te_pI       = 0;
L(n_4).Te_pI_dag   = 0;
L(n_4).Ue_pI       = 0;
L(n_4).We_pI       = 0;
L(n_4).Ve_pI       = 0;
L(n_4).Ve_pI_dag   = 0;  
L(n_4).Pvv_pI       = 0;
L(n_4).Pvv_pI_con   = 0;
L(n_4).Qvv_pI       = 0;
L(n_4).Qvv_pI_con   = 0;
L(n_4).Svv_pI       = 0;
L(n_4).Svv_pI_dag   = 0;
L(n_4).Svv_pI_con   = 0;
L(n_4).Svv_pI_dag_con = 0;
L(n_4).Rvv_pI       = 0;
L(n_4).Rvv_pI_dag   = 0;
L(n_4).Cvv_pI       = 0;
L(n_4).Cvv_pI_dag   = 0;
L(n_4).Zvv_pI       = 0;        
L(n_4).Zvv_pI_dag   = 0;
L(n_4).Qvs_pI       = 0;
L(n_4).Qvs_pI_dag   = 0;
L(n_4).Svs_pI       = 0;
L(n_4).Svs_pI_dag   = 0;
L(n_4).Svs_pI_con   = 0;
L(n_4).Svs_pI_dag_con = 0;
L(n_4).SIGvs_pI     = 0;
L(n_4).SIGvs_pI_dag = 0;
L(n_4).SIGvs_pI_con = 0;
L(n_4).SIGvs_pI_dag_con = 0;
L(n_4).Rvs_pI       = 0;
L(n_4).Rvs_pI_dag   = 0;
L(n_4).Rvs_pI_con   = 0;
L(n_4).Rvs_pI_dag_con = 0;
L(n_4).Ps_pI        = 0;
L(n_4).Ps_pI_con    = 0;
L(n_4).Cs_pI        = 0;
L(n_4).Cs_pI_dag    = 0;
%%%% Stage 3 Backward F{j=i-1} points
L(n_4).Pe_mI        = 0;
L(n_4).Pe_mI_con    = 0;
L(n_4).Ce_mI        = 0;
L(n_4).Ce_mI_dag    = 0;
L(n_4).Te_mI       = 0;
L(n_4).Te_mI_dag   = 0;
L(n_4).Ue_mI       = 0;
L(n_4).We_mI       = 0;
L(n_4).Ve_mI       = 0;
L(n_4).Ve_mI_dag   = 0;  
L(n_4).Pvv_mI       = 0;
L(n_4).Pvv_mI_con   = 0;
L(n_4).Qvv_mI       = 0;
L(n_4).Qvv_mI_con   = 0;
L(n_4).Svv_mI       = 0;
L(n_4).Svv_mI_dag   = 0;
L(n_4).Svv_mI_con   = 0;
L(n_4).Svv_mI_dag_con = 0;
L(n_4).Rvv_mI       = 0;
L(n_4).Rvv_mI_dag   = 0;
L(n_4).Cvv_mI       = 0;
L(n_4).Cvv_mI_dag   = 0;
L(n_4).Zvv_mI       = 0;        
L(n_4).Zvv_mI_dag   = 0;
L(n_4).Qvs_mI       = 0;
L(n_4).Qvs_mI_dag   = 0;
L(n_4).Svs_mI       = 0;
L(n_4).Svs_mI_dag   = 0;
L(n_4).Svs_mI_con   = 0;
L(n_4).Svs_mI_dag_con = 0;
L(n_4).SIGvs_mI     = 0;
L(n_4).SIGvs_mI_dag = 0;
L(n_4).SIGvs_mI_con = 0;
L(n_4).SIGvs_mI_dag_con = 0;
L(n_4).Rvs_mI       = 0;
L(n_4).Rvs_mI_dag   = 0;
L(n_4).Rvs_mI_con   = 0;
L(n_4).Rvs_mI_dag_con = 0;
L(n_4).Ps_mI        = 0;
L(n_4).Ps_mI_con    = 0;
L(n_4).Cs_mI        = 0;
L(n_4).Cs_mI_dag    = 0;
%%%%% Stage 4 Center F{j=i-1} Interface LEFT
L(n_4).E_ccIL        = 0;
L(n_4).E_v1IL        = 0;
L(n_4).E_v2IL        = 0;
L(n_4).E_ssIL        = 0;
L(n_4).E_cIL         = 0;
L(n_4).E_pIL         = 0;
L(n_4).E_qIL         = 0;
L(n_4).E_sIL        = 0;
L(n_4).E_rIL        = 0;
L(n_4).Pe_cIL        = 0;
L(n_4).Pe_cIL_con    = 0;
L(n_4).Ce_cIL        = 0;
L(n_4).Ce_cIL_dag    = 0;
L(n_4).Te_cIL       = 0;
L(n_4).Te_cIL_dag   = 0;
L(n_4).Ue_cIL       = 0;
L(n_4).We_cIL       = 0;
L(n_4).Ve_cIL       = 0;
L(n_4).Ve_cIL_dag   = 0;  
L(n_4).Pvv_cIL       = 0;
L(n_4).Pvv_cIL_con   = 0;
L(n_4).Qvv_cIL       = 0;
L(n_4).Qvv_cIL_con   = 0;
L(n_4).Svv_cIL       = 0;
L(n_4).Svv_cIL_dag   = 0;
L(n_4).Svv_cIL_con   = 0;
L(n_4).Svv_cIL_dag_con = 0;
L(n_4).Rvv_cIL       = 0;
L(n_4).Rvv_cIL_dag   = 0;
L(n_4).Cvv_cIL       = 0;
L(n_4).Cvv_cIL_dag   = 0;
L(n_4).Zvv_cIL       = 0;        
L(n_4).Zvv_cIL_dag   = 0;
L(n_4).Qvs_cIL       = 0;
L(n_4).Qvs_cIL_dag   = 0;
L(n_4).Svs_cIL       = 0;
L(n_4).Svs_cIL_dag   = 0;
L(n_4).Svs_cIL_con   = 0;
L(n_4).Svs_cIL_dag_con = 0;
L(n_4).SIGvs_cIL     = 0;
L(n_4).SIGvs_cIL_dag = 0;
L(n_4).SIGvs_cIL_con = 0;
L(n_4).SIGvs_cIL_dag_con = 0;
L(n_4).Rvs_cIL       = 0;
L(n_4).Rvs_cIL_dag   = 0;
L(n_4).Rvs_cIL_con   = 0;
L(n_4).Rvs_cIL_dag_con = 0;
L(n_4).Ps_cIL        = 0;
L(n_4).Ps_cIL_con    = 0;
L(n_4).Cs_cIL        = 0;
L(n_4).Cs_cIL_dag    = 0;
%%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
L(n_4).E_ccIR        = 0;
L(n_4).E_v1IR        = 0;
L(n_4).E_v2IR        = 0;
L(n_4).E_ssIR        = 0;
L(n_4).E_cIR         = 0;
L(n_4).E_pIR         = 0;
L(n_4).E_qIR         = 0;
L(n_4).E_sIR        = 0;
L(n_4).E_rIR        = 0;
L(n_4).Pe_cIR        = 0;
L(n_4).Pe_cIR_con    = 0;
L(n_4).Ce_cIR        = 0;
L(n_4).Ce_cIR_dag    = 0;
L(n_4).Te_cIR       = 0;
L(n_4).Te_cIR_dag   = 0;
L(n_4).Ue_cIR       = 0;
L(n_4).We_cIR       = 0;
L(n_4).Ve_cIR       = 0;
L(n_4).Ve_cIR_dag   = 0; 
L(n_4).Pvv_cIR       = 0;
L(n_4).Pvv_cIR_con   = 0;
L(n_4).Qvv_cIR       = 0;
L(n_4).Qvv_cIR_con   = 0;
L(n_4).Svv_cIR       = 0;
L(n_4).Svv_cIR_dag   = 0;
L(n_4).Svv_cIR_con   = 0;
L(n_4).Svv_cIR_dag_con = 0;
L(n_4).Rvv_cIR       = 0;
L(n_4).Rvv_cIR_dag   = 0;
L(n_4).Cvv_cIR       = 0;
L(n_4).Cvv_cIR_dag   = 0;
L(n_4).Zvv_cIR       = 0;        
L(n_4).Zvv_cIR_dag   = 0;
L(n_4).Qvs_cIR       = 0;
L(n_4).Qvs_cIR_dag   = 0;
L(n_4).Svs_cIR       = 0;
L(n_4).Svs_cIR_dag   = 0;
L(n_4).Svs_cIR_con   = 0;
L(n_4).Svs_cIR_dag_con = 0;
L(n_4).SIGvs_cIR     = 0;
L(n_4).SIGvs_cIR_dag = 0;
L(n_4).SIGvs_cIR_con = 0;
L(n_4).SIGvs_cIR_dag_con = 0;
L(n_4).Rvs_cIR       = 0;
L(n_4).Rvs_cIR_dag   = 0;
L(n_4).Rvs_cIR_con   = 0;
L(n_4).Rvs_cIR_dag_con = 0;
L(n_4).Ps_cIR        = 0;
L(n_4).Ps_cIR_con    = 0;
L(n_4).Cs_cIR        = 0;
L(n_4).Cs_cIR_dag    = 0;
else
%%%% Set up other interfacial terms between epilayers
L(n_4).h1 = (L(n_4).h*L(n_4-1).h*(L(n_4).h+L(n_4-1).h)); 
L(n_4).h2 = (L(n_4).h+L(n_4-1).h); 
L(n_4).h3 = (L(n_4).h*(L(n_4).h+L(n_4).h)); 
if (L(n_4).layer_cntr == 1) || (L(n_4-1).layer_cntr == 1) 
%%%% These are interface matrices of interest as they coincide with the
%%%% abrupt change in the material parameters "g".
%%%% Stage 1 Center F{j=i} points 
%%%% L(n_4).Aij_cI  = (k^2)*(((gi) + (gi-1))/2) -(-2)*(((gi)*L(n_4-1).h + (gi-1)*L(n_4).h)/L(n_4).h1); 
L(n_4).E_ccI        = ((L(n_4).E_cc + L(n_4-1).E_cc)/2);
L(n_4).E_v1I        = ((L(n_4).E_v1 + L(n_4-1).E_v1)/2);
L(n_4).E_v2I        = ((L(n_4).E_v2 + L(n_4-1).E_v2)/2);
L(n_4).E_ssI        = ((L(n_4).E_ss + L(n_4-1).E_ss)/2);
L(n_4).E_cI         = ((L(n_4).ac*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar) + L(n_4-1).ac*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar))/2);
L(n_4).E_pI         = ((L(n_4).av*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar) + L(n_4-1).av*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar))/2);
L(n_4).E_qI         = (((1/2)*L(n_4).bv*(L(n_4).euu + L(n_4).evv - 2*L(n_4).eww) + (1/2)*L(n_4-1).bv*(L(n_4-1).euu + L(n_4-1).evv - 2*L(n_4-1).eww))/2);
L(n_4).E_sI         = ((L(n_4).dv*(L(n_4).ewu - 1i*L(n_4).evw) + L(n_4-1).dv*(L(n_4-1).ewu - 1i*L(n_4-1).evw))/2);
L(n_4).E_rI         = ((((-sqrt(3)/2)*L(n_4).bv*(L(n_4).euu - L(n_4).evv) + 1i*L(n_4).dv*L(n_4).euv) + ((-sqrt(3)/2)*L(n_4-1).bv*(L(n_4-1).euu - L(n_4-1).evv) + 1i*L(n_4-1).dv*L(n_4-1).euv))/2);
L(n_4).Pe_cI        = (kpar^2)*(((L(n_4).L1eg7m) + (L(n_4-1).L1eg7m))/2)    -(-2)*(((L(n_4).L1eg7m)*L(n_4-1).h + (L(n_4-1).L1eg7m)*L(n_4).h)/L(n_4).h1);
L(n_4).Pe_cI_con    = (kpar^2)*(((L(n_4).L1eg7m) + (L(n_4-1).L1eg7m))/2)    -(-2)*(((L(n_4).L1eg7m)*L(n_4-1).h + (L(n_4-1).L1eg7m)*L(n_4).h)/L(n_4).h1);
L(n_4).Ce_cI        = 0;
L(n_4).Ce_cI_dag    = 0;
L(n_4).Te_cI       = (k_m)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/2);
L(n_4).Te_cI_dag   = (k_p)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/2);
L(n_4).Ue_cI       = 0;
L(n_4).We_cI       = 0;
L(n_4).Ve_cI       = (k_m)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/2);
L(n_4).Ve_cI_dag   = (k_p)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/2);
L(n_4).Pvv_cI       = (kpar^2)*(((L(n_4).L1vv) + (L(n_4-1).L1vv))/2)    -(-2)*(((L(n_4).L1vv)*L(n_4-1).h + (L(n_4-1).L1vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Pvv_cI_con   = (kpar^2)*(((L(n_4).L1vv) + (L(n_4-1).L1vv))/2)    -(-2)*(((L(n_4).L1vv)*L(n_4-1).h + (L(n_4-1).L1vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Qvv_cI       = (kpar^2)*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2) -2*-(-2)*(((L(n_4).L2vv)*L(n_4-1).h + (L(n_4-1).L2vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Qvv_cI_con   = (kpar^2)*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2) -2*-(-2)*(((L(n_4).L2vv)*L(n_4-1).h + (L(n_4-1).L2vv)*L(n_4).h)/L(n_4).h1);
L(n_4).Svv_cI       = 0;
L(n_4).Svv_cI_dag_con   = 0;
L(n_4).Svv_cI_con   = 0;
L(n_4).Svv_cI_dag  = 0;
L(n_4).Rvv_cI       = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vv) + (L(n_4-1).L3vv))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2); 
L(n_4).Rvv_cI_dag   = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vv) + (L(n_4-1).L3vv))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2); 
L(n_4).Cvv_cI       = 0;
L(n_4).Cvv_cI_dag   = 0;
L(n_4).Zvv_cI       = 0;        
L(n_4).Zvv_cI_dag   = 0;
L(n_4).Qvs_cI       =(kpar^2)*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2) -2*-(-2)*(((L(n_4).L2vs)*L(n_4-1).h + (L(n_4-1).L2vs)*L(n_4).h)/L(n_4).h1);
L(n_4).Qvs_cI_dag   =(kpar^2)*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2) -2*-(-2)*(((L(n_4).L2vs)*L(n_4-1).h + (L(n_4-1).L2vs)*L(n_4).h)/L(n_4).h1);
L(n_4).Svs_cI       = 0;
L(n_4).Svs_cI_dag_con   = 0;
L(n_4).Svs_cI_con       = 0;
L(n_4).Svs_cI_dag       = 0;
L(n_4).SIGvs_cI         = 0;
L(n_4).SIGvs_cI_dag_con = 0;
L(n_4).SIGvs_cI_con     = 0;
L(n_4).SIGvs_cI_dag     = 0;
L(n_4).Rvs_cI           = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_con       = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_dag_con   = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_dag       = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Ps_cI        = (kpar^2)*(((L(n_4).L1ss) + (L(n_4-1).L1ss))/2)    -(-2)*(((L(n_4).L1ss)*L(n_4-1).h + (L(n_4-1).L1ss)*L(n_4).h)/L(n_4).h1); 
L(n_4).Ps_cI_con    = (kpar^2)*(((L(n_4).L1ss) + (L(n_4-1).L1ss))/2)    -(-2)*(((L(n_4).L1ss)*L(n_4-1).h + (L(n_4-1).L1ss)*L(n_4).h)/L(n_4).h1);  
L(n_4).Cs_cI        = 0;
L(n_4).Cs_cI_dag    = 0;
%%%% Stage 2 Forward F{j=i+1) points
%%%% L(n_4).Aij_pI  =  ( )*((-2)*((((gi) +  ((gi-1)-(gi))/4)*L(n_4).h)/L(n_4).h1) + ((gi) + (gi-1))/(2*L(n_4).h2) + ((gi)/(L(n_4).h2))); 
L(n_4).Pe_pI        = ((-2)*((((L(n_4).L1eg7m) + ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Pe_pI_con    = ((-2)*((((L(n_4).L1eg7m) + ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Ce_pI        = (-1i)*(k_m)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) +((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Ce_pI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) +((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Te_pI       = 0;
L(n_4).Te_pI_dag   = 0;
L(n_4).Ue_pI       = (-1i)*(2)*(1/2)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/(2*L(n_4).h2) + ((L(n_4).Xi_cv)/(L(n_4).h2))); 
L(n_4).We_pI       = (-1i)*(1)*(1/2)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/(2*L(n_4).h2) + ((L(n_4).Xi_cs)/(L(n_4).h2))); 
L(n_4).Ve_pI       = 0;
L(n_4).Ve_pI_dag   = 0;
L(n_4).Pvv_pI       = ((-2)*((((L(n_4).L1vv) + ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Pvv_pI_con   = ((-2)*((((L(n_4).L1vv) + ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvv_pI       = (-2)*((-2)*((((L(n_4).L2vv) + ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvv_pI_con   = (-2)*((-2)*((((L(n_4).L2vv) + ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Svv_pI           = (-1i)*(sqrt(3)*k_m)*(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_pI_con       = (-1i)*(sqrt(3)*k_p)*(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Svv_pI_dag_con	= (-1i)*(sqrt(3)*k_m)*(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_pI_dag       = (-1i)*(sqrt(3)*k_p)*(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Rvv_pI       = 0; 
L(n_4).Rvv_pI_dag   = 0;
L(n_4).Cvv_pI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Cvv_pI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Zvv_pI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));        
L(n_4).Zvv_pI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Qvs_pI       = (-2)*((-2)*((((L(n_4).L2vs) + ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvs_pI_dag   = (-2)*((-2)*((((L(n_4).L2vs) + ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Svs_pI_dag_con   = (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI_dag       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI        	= (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI_con       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_dag_con = (-1i)*(sqrt(2/3)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_dag     = (-1i)*(sqrt(2/3)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI     	= (-1i)*(sqrt(2/3)*k_m)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_con     = (-1i)*(sqrt(2/3)*k_p)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Rvs_pI           = 0;
L(n_4).Rvs_pI_dag_con   = 0;
L(n_4).Rvs_pI_con       = 0;
L(n_4).Rvs_pI_dag       = 0;
L(n_4).Ps_pI        = ((-2)*((((L(n_4).L1ss) + ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Ps_pI_con    = ((-2)*((((L(n_4).L1ss) + ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Cs_pI        = (-1i)*(k_m)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss)/(L(n_4).h2)));
L(n_4).Cs_pI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss)/(L(n_4).h2)));
%%%% Stage 3 Backward F{j=i-1} points
%%%% L(n_4).Aij_pI  =   +(-2)*((((gi-1) - ((gi-1) - (gi))/4)*L(n_4-1).h)/L(n_4).h1)    - ((gi) + (gi-1))/(2*L(n_4).h2) - ((gi-1)/(L(n_4).h2)); 
L(n_4).Pe_mI        = ((-2)*((((L(n_4-1).L1eg7m) - ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Pe_mI_con    = ((-2)*((((L(n_4-1).L1eg7m) - ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Ce_mI        = (-1i)*(k_m)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Ce_mI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Te_mI       = 0;
L(n_4).Te_mI_dag   = 0;
L(n_4).Ue_mI       = (-1i)*(2)*-(1/2)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/(2*L(n_4).h2) + ((L(n_4-1).Xi_cv)/(L(n_4).h2))); 
L(n_4).We_mI       = (-1i)*(1)*-(1/2)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/(2*L(n_4).h2) + ((L(n_4-1).Xi_cs)/(L(n_4).h2))); 
L(n_4).Ve_mI       = 0;
L(n_4).Ve_mI_dag   = 0;
L(n_4).Pvv_mI       = ((-2)*((((L(n_4-1).L1vv) - ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Pvv_mI_con   = ((-2)*((((L(n_4-1).L1vv) - ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvv_mI       = (-2)*((-2)*((((L(n_4-1).L2vv) - ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvv_mI_con   = (-2)*((-2)*((((L(n_4-1).L2vv) - ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Svv_mI           = (-1i)*(sqrt(3)*k_m)*-(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_mI_con       = (-1i)*(sqrt(3)*k_p)*-(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Svv_mI_dag_con 	= (-1i)*(sqrt(3)*k_m)*-(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_mI_dag       = (-1i)*(sqrt(3)*k_p)*-(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Rvv_mI       = 0; 
L(n_4).Rvv_mI_dag   = 0;
L(n_4).Cvv_mI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Cvv_mI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Zvv_mI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));        
L(n_4).Zvv_mI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Qvs_mI       = (-2)*((-2)*((((L(n_4-1).L2vs) - ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvs_mI_dag   = (-2)*((-2)*((((L(n_4-1).L2vs) - ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Svs_mI_dag_con	= (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                            +((L(n_4-1).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI_dag       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                            +((L(n_4-1).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI        	= (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                             +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI_con       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                              +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_dag_con  = (-1i)*(sqrt(2/3)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_dag     = (-1i)*(sqrt(2/3)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI     	= (-1i)*(sqrt(2/3)*k_m)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_con     = (-1i)*(sqrt(2/3)*k_p)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                                   +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Rvs_mI           = 0;
L(n_4).Rvs_mI_dag_con   = 0;
L(n_4).Rvs_mI_con       = 0;
L(n_4).Rvs_mI_dag       = 0;
L(n_4).Ps_mI        = ((-2)*((((L(n_4-1).L1ss) - ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Ps_mI_con    = ((-2)*((((L(n_4-1).L1ss) - ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Cs_mI        = (-1i)*(k_m)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss)/(L(n_4).h2)));
L(n_4).Cs_mI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss)/(L(n_4).h2)));
%%%%% Stage 4 Center F{j=i-1} Interface LEFT 
%%%% L(n_4).Aij_cIL = (k^2)*((gi-1)/2)          -(-2)*(((3/2)*(gi-1) +(1/4)*((gi) + (gi-1)))/L(n_4).h3);  
L(n_4).E_ccIL       = (L(n_4-1).E_cc);
L(n_4).E_v1IL       = (L(n_4-1).E_v1);
L(n_4).E_v2IL       = (L(n_4-1).E_v1);
L(n_4).E_ssIL       = (L(n_4-1).E_ss);
L(n_4).E_cIL         = (L(n_4-1).ac*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar));
L(n_4).E_pIL         = (L(n_4-1).av*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar));
L(n_4).E_qIL         = ((1/2)*L(n_4-1).bv*(L(n_4-1).euu + L(n_4-1).evv - 2*L(n_4-1).eww));
L(n_4).E_sIL         = ((L(n_4-1).dv*(L(n_4-1).ewu - 1i*L(n_4-1).evw)));
L(n_4).E_rIL         = (((-sqrt(3)/2)*L(n_4-1).bv*(L(n_4-1).euu - L(n_4-1).evv) + 1i*L(n_4-1).dv*L(n_4-1).euv));
L(n_4).Pe_cIL        = (kpar^2)*(((L(n_4-1).L1eg7m)))    -(-2)*(((3/2)*(L(n_4-1).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Pe_cIL_con    = (kpar^2)*(((L(n_4-1).L1eg7m)))    -(-2)*(((3/2)*(L(n_4-1).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Ce_cIL        = 0;
L(n_4).Ce_cIL_dag    = 0;  
L(n_4).Te_cIL       = (k_m)*((L(n_4-1).Xi_cv));
L(n_4).Te_cIL_dag   = (k_p)*((L(n_4-1).Xi_cv));
L(n_4).Ue_cIL       = 0;
L(n_4).We_cIL       = 0;
L(n_4).Ve_cIL       = (k_m)*((L(n_4-1).Xi_cs));
L(n_4).Ve_cIL_dag   = (k_p)*((L(n_4-1).Xi_cs));
L(n_4).Pvv_cIL       = (kpar^2)*(((L(n_4-1).L1vv)))    -(-2)*(((3/2)*(L(n_4-1).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Pvv_cIL_con   = (kpar^2)*(((L(n_4-1).L1vv)))    -(-2)*(((3/2)*(L(n_4-1).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Qvv_cIL       = (kpar^2)*(((L(n_4-1).L2vv))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3); 
L(n_4).Qvv_cIL_con   = (kpar^2)*(((L(n_4-1).L2vv))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3);
L(n_4).Svv_cIL          = 0;
L(n_4).Svv_cIL_dag_con  = 0;
L(n_4).Svv_cIL_con      = 0;
L(n_4).Svv_cIL_dag      = 0;
L(n_4).Rvv_cIL       = (1i*2*sqrt(3)*kx*ky)*((L(n_4-1).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vv)); 
L(n_4).Rvv_cIL_dag   = (-1i*2*sqrt(3)*kx*ky)*((L(n_4-1).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vv)); 
L(n_4).Cvv_cIL       = 0;
L(n_4).Cvv_cIL_dag   = 0;
L(n_4).Zvv_cIL       = 0;        
L(n_4).Zvv_cIL_dag   = 0;
L(n_4).Qvs_cIL       = (kpar^2)*(((L(n_4-1).L2vs))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Qvs_cIL_dag   = (kpar^2)*(((L(n_4-1).L2vs))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Svs_cIL          = 0;
L(n_4).Svs_cIL_dag_con   = 0;
L(n_4).Svs_cIL_con      = 0;
L(n_4).Svs_cIL_dag      = 0;
L(n_4).SIGvs_cIL        = 0;
L(n_4).SIGvs_cIL_dag_con = 0;
L(n_4).SIGvs_cIL_con    = 0;
L(n_4).SIGvs_cIL_dag    = 0;
L(n_4).Rvs_cIL          = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_con      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_dag_con	= ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_dag      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Ps_cIL        = (kpar^2)*(((L(n_4-1).L1ss)))    -(-2)*(((3/2)*(L(n_4-1).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3); 
L(n_4).Ps_cIL_con    = (kpar^2)*(((L(n_4-1).L1ss)))    -(-2)*(((3/2)*(L(n_4-1).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3);  
L(n_4).Cs_cIL        = 0;
L(n_4).Cs_cIL_dag    = 0;    
%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
%%%% L(n_4).Aij_cIR = (k^2)*((gi)/2)            -(-2)*(((3/2)*(gi) + (1/4)*((gi) + (gi-1)))/L(n_4).h3); 
L(n_4).E_ccIR = (L(n_4).E_cc);
L(n_4).E_v1IR = (L(n_4).E_v1);
L(n_4).E_v2IR = (L(n_4).E_v2);
L(n_4).E_ssIR = (L(n_4).E_ss); 
L(n_4).E_cIR         = (L(n_4).ac*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_pIR         = (L(n_4).av*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_qIR         = ((1/2)*L(n_4).bv*(L(n_4).euu + L(n_4).evv - 2*L(n_4).eww));
L(n_4).E_sIR         = (L(n_4).dv*(L(n_4).ewu - 1i*L(n_4).evw));
L(n_4).E_rIR         = (((-sqrt(3)/2)*L(n_4).bv*(L(n_4).euu - L(n_4).evv) + 1i*L(n_4).dv*L(n_4).euv));
L(n_4).Pe_cIR        = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((3/2)*(L(n_4).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Pe_cIR_con    = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((3/2)*(L(n_4).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Ce_cIR        = 0;
L(n_4).Ce_cIR_dag    = 0; 
L(n_4).Te_cIR       = (k_m)*((L(n_4).Xi_cv));
L(n_4).Te_cIR_dag   = (k_p)*((L(n_4).Xi_cv));
L(n_4).Ue_cIR       = 0;
L(n_4).We_cIR       = 0;
L(n_4).Ve_cIR       = (k_m)*((L(n_4).Xi_cs));
L(n_4).Ve_cIR_dag   = (k_p)*((L(n_4).Xi_cs));
L(n_4).Pvv_cIR       = (kpar^2)*(((L(n_4).L1vv)))    -(-2)*(((3/2)*(L(n_4).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Pvv_cIR_con   = (kpar^2)*(((L(n_4).L1vv)))    -(-2)*(((3/2)*(L(n_4).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Qvv_cIR       = (kpar^2)*(((L(n_4).L2vv))) -2*-(-2)*(((3/2)*(L(n_4).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3); 
L(n_4).Qvv_cIR_con   = (kpar^2)*(((L(n_4).L2vv))) -2*-(-2)*(((3/2)*(L(n_4).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3);
L(n_4).Svv_cIR          = 0;
L(n_4).Svv_cIR_dag_con 	= 0;
L(n_4).Svv_cIR_con      = 0;
L(n_4).Svv_cIR_dag      = 0;
L(n_4).Rvv_cIR       = ( 1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Rvv_cIR_dag   = (-1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Cvv_cIR       = 0;
L(n_4).Cvv_cIR_dag   = 0;
L(n_4).Zvv_cIR       = 0;        
L(n_4).Zvv_cIR_dag   = 0;
L(n_4).Qvs_cIR       = (kpar^2)*(((L(n_4).L2vs))) -2*-(-2)*(((3/2)*(L(n_4).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Qvs_cIR_dag   = (kpar^2)*(((L(n_4).L2vs))) -2*-(-2)*(((3/2)*(L(n_4).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Svs_cIR          = 0;
L(n_4).Svs_cIR_dag_con   = 0;
L(n_4).Svs_cIR_con      = 0;
L(n_4).Svs_cIR_dag      = 0;
L(n_4).SIGvs_cIR        = 0;
L(n_4).SIGvs_cIR_dag_con = 0;
L(n_4).SIGvs_cIR_con    = 0;
L(n_4).SIGvs_cIR_dag    = 0;
L(n_4).Rvs_cIR          = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_con      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_dag_con 	= ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_dag      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Ps_cIR        = (kpar^2)*(((L(n_4).L1ss)))    -(-2)*(((3/2)*(L(n_4).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3); 
L(n_4).Ps_cIR_con    = (kpar^2)*(((L(n_4).L1ss)))    -(-2)*(((3/2)*(L(n_4).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3);  
L(n_4).Cs_cIR        = 0;
L(n_4).Cs_cIR_dag    = 0;         
else
%%%% Define all other quasi-interfacial points within the heterostructure
%%%% Stage 1 Center F{j=i} points 
%%%% L(n_4).Aij_cI  = (k^2)*(((gi) + (gi-1))/2) -(-2)*(((gi)*L(n_4-1).h + (gi-1)*L(n_4).h)/L(n_4).h1); 
L(n_4).E_ccI        = ((L(n_4).E_cc + L(n_4-1).E_cc)/2);
L(n_4).E_v1I        = ((L(n_4).E_v1 + L(n_4-1).E_v1)/2);
L(n_4).E_v2I        = ((L(n_4).E_v2 + L(n_4-1).E_v2)/2);
L(n_4).E_ssI        = ((L(n_4).E_ss + L(n_4-1).E_ss)/2);
L(n_4).E_cI         = ((L(n_4).ac*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar) + L(n_4-1).ac*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar))/2);
L(n_4).E_pI         = ((L(n_4).av*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar) + L(n_4-1).av*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar))/2);
L(n_4).E_qI         = (((1/2)*L(n_4).bv*(L(n_4).euu + L(n_4).evv - 2*L(n_4).eww) + (1/2)*L(n_4-1).bv*(L(n_4-1).euu + L(n_4-1).evv - 2*L(n_4-1).eww))/2);
L(n_4).E_sI         = ((L(n_4).dv*(L(n_4).ewu - 1i*L(n_4).evw) + L(n_4-1).dv*(L(n_4-1).ewu - 1i*L(n_4-1).evw))/2);
L(n_4).E_rI         = ((((-sqrt(3)/2)*L(n_4).bv*(L(n_4).euu - L(n_4).evv) + 1i*L(n_4).dv*L(n_4).euv) + ((-sqrt(3)/2)*L(n_4-1).bv*(L(n_4-1).euu - L(n_4-1).evv) + 1i*L(n_4-1).dv*L(n_4-1).euv))/2);
L(n_4).Pe_cI        = (kpar^2)*(((L(n_4).L1eg7m) + (L(n_4-1).L1eg7m))/2)    -(-2)*(((L(n_4).L1eg7m)*L(n_4-1).h + (L(n_4-1).L1eg7m)*L(n_4).h)/L(n_4).h1);
L(n_4).Pe_cI_con    = (kpar^2)*(((L(n_4).L1eg7m) + (L(n_4-1).L1eg7m))/2)    -(-2)*(((L(n_4).L1eg7m)*L(n_4-1).h + (L(n_4-1).L1eg7m)*L(n_4).h)/L(n_4).h1);
L(n_4).Ce_cI        = 0;
L(n_4).Ce_cI_dag    = 0;
L(n_4).Te_cI       = (k_m)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/2);
L(n_4).Te_cI_dag   = (k_p)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/2);
L(n_4).Ue_cI       = 0;
L(n_4).We_cI       = 0;
L(n_4).Ve_cI       = (k_m)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/2);
L(n_4).Ve_cI_dag   = (k_p)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/2);
L(n_4).Pvv_cI       = (kpar^2)*(((L(n_4).L1vv) + (L(n_4-1).L1vv))/2)    -(-2)*(((L(n_4).L1vv)*L(n_4-1).h + (L(n_4-1).L1vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Pvv_cI_con   = (kpar^2)*(((L(n_4).L1vv) + (L(n_4-1).L1vv))/2)    -(-2)*(((L(n_4).L1vv)*L(n_4-1).h + (L(n_4-1).L1vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Qvv_cI       = (kpar^2)*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2) -2*-(-2)*(((L(n_4).L2vv)*L(n_4-1).h + (L(n_4-1).L2vv)*L(n_4).h)/L(n_4).h1); 
L(n_4).Qvv_cI_con   = (kpar^2)*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2) -2*-(-2)*(((L(n_4).L2vv)*L(n_4-1).h + (L(n_4-1).L2vv)*L(n_4).h)/L(n_4).h1);
L(n_4).Svv_cI       = 0;
L(n_4).Svv_cI_dag_con   = 0;
L(n_4).Svv_cI_con   = 0;
L(n_4).Svv_cI_dag  = 0;
L(n_4).Rvv_cI       = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vv) + (L(n_4-1).L3vv))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2); 
L(n_4).Rvv_cI_dag   = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vv) + (L(n_4-1).L3vv))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vv) + (L(n_4-1).L2vv))/2); 
L(n_4).Cvv_cI       = 0;
L(n_4).Cvv_cI_dag   = 0;
L(n_4).Zvv_cI       = 0;        
L(n_4).Zvv_cI_dag   = 0;
L(n_4).Qvs_cI       =(kpar^2)*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2) -2*-(-2)*(((L(n_4).L2vs)*L(n_4-1).h + (L(n_4-1).L2vs)*L(n_4).h)/L(n_4).h1);
L(n_4).Qvs_cI_dag   =(kpar^2)*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2) -2*-(-2)*(((L(n_4).L2vs)*L(n_4-1).h + (L(n_4-1).L2vs)*L(n_4).h)/L(n_4).h1);
L(n_4).Svs_cI       = 0;
L(n_4).Svs_cI_dag_con   = 0;
L(n_4).Svs_cI_con       = 0;
L(n_4).Svs_cI_dag       = 0;
L(n_4).SIGvs_cI         = 0;
L(n_4).SIGvs_cI_dag_con = 0;
L(n_4).SIGvs_cI_con     = 0;
L(n_4).SIGvs_cI_dag     = 0;
L(n_4).Rvs_cI           = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_con       = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_dag_con   = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Rvs_cI_dag       = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs) + (L(n_4-1).L3vs))/2) - (sqrt(3)*(kx*kx - ky*ky))*(((L(n_4).L2vs) + (L(n_4-1).L2vs))/2);
L(n_4).Ps_cI        = (kpar^2)*(((L(n_4).L1ss) + (L(n_4-1).L1ss))/2)    -(-2)*(((L(n_4).L1ss)*L(n_4-1).h + (L(n_4-1).L1ss)*L(n_4).h)/L(n_4).h1); 
L(n_4).Ps_cI_con    = (kpar^2)*(((L(n_4).L1ss) + (L(n_4-1).L1ss))/2)    -(-2)*(((L(n_4).L1ss)*L(n_4-1).h + (L(n_4-1).L1ss)*L(n_4).h)/L(n_4).h1);  
L(n_4).Cs_cI        = 0;
L(n_4).Cs_cI_dag    = 0;
%%%% Stage 2 Forward F{j=i+1) points
%%%% L(n_4).Aij_pI  =  ( )*((-2)*((((gi) +  ((gi-1)-(gi))/4)*L(n_4).h)/L(n_4).h1) + ((gi) + (gi-1))/(2*L(n_4).h2) + ((gi)/(L(n_4).h2))); 
L(n_4).Pe_pI        = ((-2)*((((L(n_4).L1eg7m) + ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Pe_pI_con    = ((-2)*((((L(n_4).L1eg7m) + ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Ce_pI        = (-1i)*(k_m)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Ce_pI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Te_pI       = 0;
L(n_4).Te_pI_dag   = 0;
L(n_4).Ue_pI       = (-1i)*(2)*(1/2)*((((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/(2*L(n_4).h2)) + ((L(n_4).Xi_cv)/(L(n_4).h2))); 
L(n_4).We_pI       = (-1i)*(1)*(1/2)*((((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/(2*L(n_4).h2)) + ((L(n_4).Xi_cs)/(L(n_4).h2))); 
L(n_4).Ve_pI       = 0;
L(n_4).Ve_pI_dag   = 0;
L(n_4).Pvv_pI       = ((-2)*((((L(n_4).L1vv) + ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Pvv_pI_con   = ((-2)*((((L(n_4).L1vv) + ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvv_pI       = (-2)*((-2)*((((L(n_4).L2vv) + ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvv_pI_con   = (-2)*((-2)*((((L(n_4).L2vv) + ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Svv_pI           = (-1i)*(sqrt(3)*k_m)*(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_pI_con       = (-1i)*(sqrt(3)*k_p)*(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Svv_pI_dag_con	= (-1i)*(sqrt(3)*k_m)*(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_pI_dag       = (-1i)*(sqrt(3)*k_p)*(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Rvv_pI       = 0; 
L(n_4).Rvv_pI_dag   = 0;
L(n_4).Cvv_pI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Cvv_pI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Zvv_pI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv)/(L(n_4).h2)));        
L(n_4).Zvv_pI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Qvs_pI       = (-2)*((-2)*((((L(n_4).L2vs) + ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Qvs_pI_dag   = (-2)*((-2)*((((L(n_4).L2vs) + ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Svs_pI_dag_con 	= (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI_dag       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI        	= (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Svs_pI_con       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_dag_con = (-1i)*(sqrt(2/3)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_dag     = (-1i)*(sqrt(2/3)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI     	= (-1i)*(sqrt(2/3)*k_m)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_pI_con     = (-1i)*(sqrt(2/3)*k_p)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Rvs_pI           = 0;
L(n_4).Rvs_pI_dag_con   = 0;
L(n_4).Rvs_pI_con       = 0;
L(n_4).Rvs_pI_dag       = 0;
L(n_4).Ps_pI        = ((-2)*((((L(n_4).L1ss) + ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4).h)/L(n_4).h1)); 
L(n_4).Ps_pI_con    = ((-2)*((((L(n_4).L1ss) + ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4).h)/L(n_4).h1));
L(n_4).Cs_pI        = (-1i)*(k_m)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss)/(L(n_4).h2)));
L(n_4).Cs_pI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss)/(L(n_4).h2)));
%%%% Stage 3 Backward F{j=i-1} points
%%%% L(n_4).Aij_pI  =   +(-2)*((((gi-1) - ((gi-1) - (gi))/4)*L(n_4-1).h)/L(n_4).h1)    - ((gi) + (gi-1))/(2*L(n_4).h2) - ((gi-1)/(L(n_4).h2)); 
L(n_4).Pe_mI        = ((-2)*((((L(n_4-1).L1eg7m) - ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Pe_mI_con    = ((-2)*((((L(n_4-1).L1eg7m) - ((L(n_4-1).L1eg7m) - (L(n_4).L1eg7m))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Ce_mI        = (-1i)*(k_m)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Ce_mI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7p_cc - 2*L(n_4).Z_g8p_cc) + (L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7p_cc - 2*L(n_4-1).Z_g8p_cc)/(L(n_4).h2)));
L(n_4).Te_mI       = 0;
L(n_4).Te_mI_dag   = 0;
L(n_4).Ue_mI       = (-1i)*(2)*-(1/2)*(((L(n_4).Xi_cv) + (L(n_4-1).Xi_cv))/(2*L(n_4).h2) + ((L(n_4-1).Xi_cv)/(L(n_4).h2))); 
L(n_4).We_mI       = (-1i)*(1)*-(1/2)*(((L(n_4).Xi_cs) + (L(n_4-1).Xi_cs))/(2*L(n_4).h2) + ((L(n_4-1).Xi_cs)/(L(n_4).h2))); 
L(n_4).Ve_mI       = 0;
L(n_4).Ve_mI_dag   = 0;
L(n_4).Pvv_mI       = ((-2)*((((L(n_4-1).L1vv) - ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Pvv_mI_con   = ((-2)*((((L(n_4-1).L1vv) - ((L(n_4-1).L1vv) - (L(n_4).L1vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvv_mI       = (-2)*((-2)*((((L(n_4-1).L2vv) - ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvv_mI_con   = (-2)*((-2)*((((L(n_4-1).L2vv) - ((L(n_4-1).L2vv) - (L(n_4).L2vv))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Svv_mI           = (-1i)*(sqrt(3)*k_m)*-(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_mI_con       = (-1i)*(sqrt(3)*k_p)*-(((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Svv_mI_dag_con 	= (-1i)*(sqrt(3)*k_m)*-(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                             +((2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Svv_mI_dag       = (-1i)*(sqrt(3)*k_p)*-(((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - 4*L(n_4-1).Z_g8m2_vv - L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g6m_vv + L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Rvv_mI       = 0; 
L(n_4).Rvv_mI_dag   = 0;
L(n_4).Cvv_mI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Cvv_mI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g7m_vv - L(n_4).Z_g8m1_vv - 4*L(n_4).Z_g8m2_vv - 5*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g7m_vv - L(n_4-1).Z_g8m1_vv - 4*L(n_4-1).Z_g8m2_vv - 5*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));
L(n_4).Zvv_mI       = (-1i)*(k_m)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2)));        
L(n_4).Zvv_mI_dag   = (-1i)*(k_p)*(-((2*L(n_4).Z_g6m_vv - L(n_4).Z_g8m1_vv - 3*L(n_4).Z_g8m3_vv) + (2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv))/(2*L(n_4).h2) ...
                     +((2*L(n_4-1).Z_g6m_vv - L(n_4-1).Z_g8m1_vv - 3*L(n_4-1).Z_g8m3_vv)/(L(n_4).h2))); 
L(n_4).Qvs_mI       = (-2)*((-2)*((((L(n_4-1).L2vs) - ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Qvs_mI_dag   = (-2)*((-2)*((((L(n_4-1).L2vs) - ((L(n_4-1).L2vs) - (L(n_4).L2vs))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Svs_mI_dag_con  = (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                            +((L(n_4-1).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI_dag       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                            +((L(n_4-1).Z_g8m1_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI        	= (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                             +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Svs_mI_con       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g8m1_vs) + (L(n_4-1).Z_g8m1_vs))/(2*L(n_4).h2) ...
                              +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_dag_con = (-1i)*(sqrt(2/3)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_dag     = (-1i)*(sqrt(2/3)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs) + (L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI     	= (-1i)*(sqrt(2/3)*k_m)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                               +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).SIGvs_mI_con     = (-1i)*(sqrt(2/3)*k_p)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs) + (2*L(n_4-1).Z_g7m_vs + L(n_4-1).Z_g8m1_vs + 4*L(n_4-1).Z_g8m2_vs))/(2*L(n_4).h2) ...
                                   +((L(n_4-1).Z_g7m_vs + 2*L(n_4-1).Z_g8m1_vs + 2*L(n_4-1).Z_g8m2_vs)/(L(n_4).h2))); 
L(n_4).Rvs_mI           = 0;
L(n_4).Rvs_mI_dag_con   = 0;
L(n_4).Rvs_mI_con       = 0;
L(n_4).Rvs_mI_dag       = 0;
L(n_4).Ps_mI        = ((-2)*((((L(n_4-1).L1ss) - ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4-1).h)/L(n_4).h1)); 
L(n_4).Ps_mI_con    = ((-2)*((((L(n_4-1).L1ss) - ((L(n_4-1).L1ss) - (L(n_4).L1ss))/4)*L(n_4-1).h)/L(n_4).h1));
L(n_4).Cs_mI        = (-1i)*(k_m)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss)/(L(n_4).h2)));
L(n_4).Cs_mI_dag    = (-1i)*(k_p)*(-((L(n_4).Z_g7m_ss - 2*L(n_4).Z_g8m_ss) + (L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss))/(2*L(n_4).h2) ...
                     +((L(n_4-1).Z_g7m_ss - 2*L(n_4-1).Z_g8m_ss)/(L(n_4).h2)));
%%%% Stage 4 Center F{j=i-1} Interface LEFT 
%%%% L(n_4).Aij_cIL = (k^2)*((gi-1)/2)          -(-2)*(((3/2)*(gi-1) + (1/4)*((gi) + (gi-1)))/L(n_4).h3);  
L(n_4).E_ccIL       = (L(n_4-1).E_cc);
L(n_4).E_v1IL       = (L(n_4-1).E_v1);
L(n_4).E_v2IL       = (L(n_4-1).E_v1);
L(n_4).E_ssIL       = (L(n_4-1).E_ss);
L(n_4).E_cIL         = (L(n_4-1).ac*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar));
L(n_4).E_pIL         = (L(n_4-1).av*(L(n_4-1).euu + L(n_4-1).evv + L(n_4-1).eww - 3*L(n_4-1).epar));
L(n_4).E_qIL         = ((1/2)*L(n_4-1).bv*(L(n_4-1).euu + L(n_4-1).evv - 2*L(n_4-1).eww));
L(n_4).E_sIL         = ((L(n_4-1).dv*(L(n_4-1).ewu - 1i*L(n_4-1).evw)));
L(n_4).E_rIL         = (((-sqrt(3)/2)*L(n_4-1).bv*(L(n_4-1).euu - L(n_4-1).evv) + 1i*L(n_4-1).dv*L(n_4-1).euv));
L(n_4).Pe_cIL        = (kpar^2)*(((L(n_4-1).L1eg7m)))    -(-2)*(((3/2)*(L(n_4-1).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Pe_cIL_con    = (kpar^2)*(((L(n_4-1).L1eg7m)))    -(-2)*(((3/2)*(L(n_4-1).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Ce_cIL        = 0;
L(n_4).Ce_cIL_dag    = 0;
L(n_4).Te_cIL       = (k_m)*((L(n_4-1).Xi_cv));
L(n_4).Te_cIL_dag   = (k_p)*((L(n_4-1).Xi_cv));
L(n_4).Ue_cIL       = 0;
L(n_4).We_cIL       = 0;
L(n_4).Ve_cIL       = (k_m)*((L(n_4-1).Xi_cs));
L(n_4).Ve_cIL_dag   = (k_p)*((L(n_4-1).Xi_cs));
L(n_4).Pvv_cIL       = (kpar^2)*(((L(n_4-1).L1vv)))    -(-2)*(((3/2)*(L(n_4-1).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Pvv_cIL_con   = (kpar^2)*(((L(n_4-1).L1vv)))    -(-2)*(((3/2)*(L(n_4-1).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Qvv_cIL       = (kpar^2)*(((L(n_4-1).L2vv))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3); 
L(n_4).Qvv_cIL_con   = (kpar^2)*(((L(n_4-1).L2vv))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3);
L(n_4).Svv_cIL          = 0;
L(n_4).Svv_cIL_dag_con  = 0;
L(n_4).Svv_cIL_con      = 0;
L(n_4).Svv_cIL_dag      = 0;
L(n_4).Rvv_cIL       = (1i*2*sqrt(3)*kx*ky)*((L(n_4-1).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vv)); 
L(n_4).Rvv_cIL_dag   = (-1i*2*sqrt(3)*kx*ky)*((L(n_4-1).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vv)); 
L(n_4).Cvv_cIL       = 0;
L(n_4).Cvv_cIL_dag   = 0;
L(n_4).Zvv_cIL       = 0;        
L(n_4).Zvv_cIL_dag   = 0;
L(n_4).Qvs_cIL       = (kpar^2)*(((L(n_4-1).L2vs))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Qvs_cIL_dag   = (kpar^2)*(((L(n_4-1).L2vs))) -2*-(-2)*(((3/2)*(L(n_4-1).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Svs_cIL          = 0;
L(n_4).Svs_cIL_dag_con   = 0;
L(n_4).Svs_cIL_con      = 0;
L(n_4).Svs_cIL_dag      = 0;
L(n_4).SIGvs_cIL        = 0;
L(n_4).SIGvs_cIL_dag_con = 0;
L(n_4).SIGvs_cIL_con    = 0;
L(n_4).SIGvs_cIL_dag    = 0;
L(n_4).Rvs_cIL          = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_con      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_dag_con	= ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Rvs_cIL_dag      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4-1).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4-1).L2vs));
L(n_4).Ps_cIL        = (kpar^2)*(((L(n_4-1).L1ss)))    -(-2)*(((3/2)*(L(n_4-1).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3); 
L(n_4).Ps_cIL_con    = (kpar^2)*(((L(n_4-1).L1ss)))    -(-2)*(((3/2)*(L(n_4-1).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3);  
L(n_4).Cs_cIL        = 0;
L(n_4).Cs_cIL_dag    = 0;    
%%%% Stage 5 Center F{j=i+1} Interface RIGHT 
%%%% L(n_4).Aij_cIR = (k^2)*((gi)/2)            -(-2)*(((3/2)*(gi) + (1/4)*((gi) + (gi-1)))/L(n_4).h3); 
L(n_4).E_ccIR = (L(n_4).E_cc);
L(n_4).E_v1IR = (L(n_4).E_v1);
L(n_4).E_v2IR = (L(n_4).E_v2);
L(n_4).E_ssIR = (L(n_4).E_ss); 
L(n_4).E_cIR         = (L(n_4).ac*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_pIR         = (L(n_4).av*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_qIR         = ((1/2)*L(n_4).bv*(L(n_4).euu + L(n_4).evv - 2*L(n_4).eww));
L(n_4).E_sIR         = (L(n_4).dv*(L(n_4).ewu - 1i*L(n_4).evw));
L(n_4).E_rIR         = (((-sqrt(3)/2)*L(n_4).bv*(L(n_4).euu - L(n_4).evv) + 1i*L(n_4).dv*L(n_4).euv));
L(n_4).Pe_cIR        = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((3/2)*(L(n_4).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Pe_cIR_con    = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((3/2)*(L(n_4).L1eg7m) + (1/4)*(L(n_4).L1eg7m + L(n_4-1).L1eg7m))/L(n_4).h3);
L(n_4).Ce_cIR        = 0;
L(n_4).Ce_cIR_dag    = 0;  
L(n_4).Te_cIR       = (k_m)*((L(n_4).Xi_cv));
L(n_4).Te_cIR_dag   = (k_p)*((L(n_4).Xi_cv));
L(n_4).Ue_cIR       = 0;
L(n_4).We_cIR       = 0;
L(n_4).Ve_cIR       = (k_m)*((L(n_4).Xi_cs));
L(n_4).Ve_cIR_dag   = (k_p)*((L(n_4).Xi_cs));
L(n_4).Pvv_cIR       = (kpar^2)*(((L(n_4).L1vv)))    -(-2)*(((3/2)*(L(n_4).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Pvv_cIR_con   = (kpar^2)*(((L(n_4).L1vv)))    -(-2)*(((3/2)*(L(n_4).L1vv) + (1/4)*(L(n_4).L1vv + L(n_4-1).L1vv))/L(n_4).h3); 
L(n_4).Qvv_cIR       = (kpar^2)*(((L(n_4).L2vv))) -2*-(-2)*(((3/2)*(L(n_4).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3); 
L(n_4).Qvv_cIR_con   = (kpar^2)*(((L(n_4).L2vv))) -2*-(-2)*(((3/2)*(L(n_4).L2vv) + (1/4)*(L(n_4).L2vv + L(n_4-1).L2vv))/L(n_4).h3);
L(n_4).Svv_cIR          = 0;
L(n_4).Svv_cIR_dag_con 	= 0;
L(n_4).Svv_cIR_con      = 0;
L(n_4).Svv_cIR_dag      = 0;
L(n_4).Rvv_cIR       = ( 1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Rvv_cIR_dag   = (-1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Cvv_cIR       = 0;
L(n_4).Cvv_cIR_dag   = 0;
L(n_4).Zvv_cIR       = 0;        
L(n_4).Zvv_cIR_dag   = 0;
L(n_4).Qvs_cIR       = (kpar^2)*(((L(n_4).L2vs))) -2*-(-2)*(((3/2)*(L(n_4).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Qvs_cIR_dag   = (kpar^2)*(((L(n_4).L2vs))) -2*-(-2)*(((3/2)*(L(n_4).L2vs) + (1/4)*(L(n_4).L2vs + L(n_4-1).L2vs))/L(n_4).h3);
L(n_4).Svs_cIR          = 0;
L(n_4).Svs_cIR_dag_con   = 0;
L(n_4).Svs_cIR_con      = 0;
L(n_4).Svs_cIR_dag      = 0;
L(n_4).SIGvs_cIR        = 0;
L(n_4).SIGvs_cIR_dag_con = 0;
L(n_4).SIGvs_cIR_con    = 0;
L(n_4).SIGvs_cIR_dag    = 0;
L(n_4).Rvs_cIR          = ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_con      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_dag_con 	= ( 1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Rvs_cIR_dag      = (-1i*2*sqrt(3)*kx*ky)*(((L(n_4).L3vs))) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vs));
L(n_4).Ps_cIR        = (kpar^2)*(((L(n_4).L1ss)))    -(-2)*(((3/2)*(L(n_4).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3); 
L(n_4).Ps_cIR_con    = (kpar^2)*(((L(n_4).L1ss)))    -(-2)*(((3/2)*(L(n_4).L1ss) + (1/4)*(L(n_4).L1ss + L(n_4-1).L1ss))/L(n_4).h3);  
L(n_4).Cs_cIR        = 0;
L(n_4).Cs_cIR_dag    = 0;         
end
end
%%%% Set up terms in the epilayer away from the interfacial mesh points 
L(n_4).h4 = ((L(n_4).h*L(n_4).h*(L(n_4).h+L(n_4).h)));
L(n_4).h5 = ((L(n_4).h + L(n_4).h)); 
%%%% Stage 1 Center F{j=i} points 
%%%% L(n_4).Aij_cI  = (k^2)*(((gi) + (gi-1))/2) -(-2)*(((gi)*L(n_4-1).h + (gi-1)*L(n_4).h)/L(n_4).h1); 
L(n_4).E_ccJ = (L(n_4).E_cc);
L(n_4).E_v1J = (L(n_4).E_v1);
L(n_4).E_v2J = (L(n_4).E_v2);
L(n_4).E_ssJ = (L(n_4).E_ss);
L(n_4).E_cJ  = (L(n_4).ac*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_pJ  = (L(n_4).av*(L(n_4).euu + L(n_4).evv + L(n_4).eww - 3*L(n_4).epar));
L(n_4).E_qJ  = ((1/2)*L(n_4).bv*(L(n_4).euu + L(n_4).evv - 2*L(n_4).eww));
L(n_4).E_sJ  = (L(n_4).dv*(L(n_4).ewu - 1i*L(n_4).evw));
L(n_4).E_rJ  = (((-sqrt(3)/2)*L(n_4).bv*(L(n_4).euu - L(n_4).evv) + 1i*L(n_4).dv*L(n_4).euv));
L(n_4).Pe_cJ        = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((L(n_4).L1eg7m)*L(n_4).h + (L(n_4).L1eg7m)*L(n_4).h)/L(n_4).h4);
L(n_4).Pe_cJ_con    = (kpar^2)*(((L(n_4).L1eg7m)))    -(-2)*(((L(n_4).L1eg7m)*L(n_4).h + (L(n_4).L1eg7m)*L(n_4).h)/L(n_4).h4);
L(n_4).Ce_cJ        = 0;
L(n_4).Ce_cJ_dag    = 0;
L(n_4).Te_cJ       = (k_m)*((L(n_4).Xi_cv));
L(n_4).Te_cJ_dag   = (k_p)*((L(n_4).Xi_cv));
L(n_4).Ue_cJ       = 0;
L(n_4).We_cJ       = 0;
L(n_4).Ve_cJ       = (k_m)*((L(n_4).Xi_cs));
L(n_4).Ve_cJ_dag   = (k_p)*((L(n_4).Xi_cs));
L(n_4).Pvv_cJ       = (kpar^2)*(L(n_4).L1vv)    -(-2)*(((L(n_4).L1vv)*L(n_4).h + (L(n_4).L1vv)*L(n_4).h)/L(n_4).h4); 
L(n_4).Pvv_cJ_con   = (kpar^2)*(L(n_4).L1vv)    -(-2)*(((L(n_4).L1vv)*L(n_4).h + (L(n_4).L1vv)*L(n_4).h)/L(n_4).h4); 
L(n_4).Qvv_cJ       = (kpar^2)*(L(n_4).L2vv) -2*-(-2)*(((L(n_4).L2vv)*L(n_4).h + (L(n_4).L2vv)*L(n_4).h)/L(n_4).h4); 
L(n_4).Qvv_cJ_con   = (kpar^2)*(L(n_4).L2vv) -2*-(-2)*(((L(n_4).L2vv)*L(n_4).h + (L(n_4).L2vv)*L(n_4).h)/L(n_4).h4);
L(n_4).Svv_cJ           = 0;
L(n_4).Svv_cJ_con       = 0;
L(n_4).Svv_cJ_dag_con 	= 0;
L(n_4).Svv_cJ_dag       = 0;
L(n_4).Rvv_cJ       = (1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Rvv_cJ_dag   = (-1i*2*sqrt(3)*kx*ky)*((L(n_4).L3vv)) - (sqrt(3)*(kx*kx - ky*ky))*((L(n_4).L2vv)); 
L(n_4).Cvv_cJ       = 0;
L(n_4).Cvv_cJ_dag   = 0;
L(n_4).Zvv_cJ       = 0;        
L(n_4).Zvv_cJ_dag   = 0;
L(n_4).Qvs_cJ       = (kpar^2)*(L(n_4).L2vs) -2*-(-2)*(((L(n_4).L2vs)*L(n_4).h + (L(n_4).L2vs)*L(n_4).h)/L(n_4).h4);
L(n_4).Qvs_cJ_dag   = (kpar^2)*(L(n_4).L2vs) -2*-(-2)*(((L(n_4).L2vs)*L(n_4).h + (L(n_4).L2vs)*L(n_4).h)/L(n_4).h4);
L(n_4).Svs_cJ           = 0;
L(n_4).Svs_cJ_con       = 0;
L(n_4).Svs_cJ_dag_con 	= 0;
L(n_4).Svs_cJ_dag       = 0;
L(n_4).SIGvs_cJ         = 0;
L(n_4).SIGvs_cJ_con     = 0;
L(n_4).SIGvs_cJ_dag_con = 0;
L(n_4).SIGvs_cJ_dag     = 0;
L(n_4).Rvs_cJ           = ( 1i*2*sqrt(3)*kx*ky)*(L(n_4).L3vs) - (sqrt(3)*(kx*kx - ky*ky))*(L(n_4).L2vs);
L(n_4).Rvs_cJ_con       = (-1i*2*sqrt(3)*kx*ky)*(L(n_4).L3vs) - (sqrt(3)*(kx*kx - ky*ky))*(L(n_4).L2vs);
L(n_4).Rvs_cJ_dag_con   = ( 1i*2*sqrt(3)*kx*ky)*(L(n_4).L3vs) - (sqrt(3)*(kx*kx - ky*ky))*(L(n_4).L2vs);
L(n_4).Rvs_cJ_dag       = (-1i*2*sqrt(3)*kx*ky)*(L(n_4).L3vs) - (sqrt(3)*(kx*kx - ky*ky))*(L(n_4).L2vs);
L(n_4).Ps_cJ      	= (kpar^2)*(L(n_4).L1ss)    -(-2)*(((L(n_4).L1ss)*L(n_4).h + (L(n_4).L1ss)*L(n_4).h)/L(n_4).h4); 
L(n_4).Ps_cJ_con    = (kpar^2)*(L(n_4).L1ss)    -(-2)*(((L(n_4).L1ss)*L(n_4).h + (L(n_4).L1ss)*L(n_4).h)/L(n_4).h4);  
L(n_4).Cs_cJ        = 0;
L(n_4).Cs_cJ_dag    = 0;
%%%% Stage 2 Forward F{j=i+1) points
%%%% L(n_4).Aij_pI  =  ( )*((-2)*((((gi) +  ((gi-1) -(gi))/4)*L(n_4).h)/L(n_4).h1) + ((gi) + (gi-1))/(2*L(n_4).h2) + ((gi)/(L(n_4).h2))); 
L(n_4).Pe_pJ        = ((-2)*((((L(n_4).L1eg7m))*L(n_4).h)/L(n_4).h4)); 
L(n_4).Pe_pJ_con    = ((-2)*((((L(n_4).L1eg7m))*L(n_4).h)/L(n_4).h4)); 
L(n_4).Ce_pJ        = 0;
L(n_4).Ce_pJ_dag    = 0;
L(n_4).Te_pJ       = 0;
L(n_4).Te_pJ_dag   = 0;
L(n_4).Ue_pJ       = (-1i)*(2)*(1/2)*(((L(n_4).Xi_cv)/(L(n_4).h5)) + ((L(n_4).Xi_cv)/(L(n_4).h5))); 
L(n_4).We_pJ       = (-1i)*(1)*(1/2)*(((L(n_4).Xi_cs)/(L(n_4).h5)) + ((L(n_4).Xi_cs)/(L(n_4).h5))); 
L(n_4).Ve_pJ       = 0;
L(n_4).Ve_pJ_dag   = 0;
L(n_4).Pvv_pJ       = ((-2)*(((L(n_4).L1vv)*L(n_4).h)/L(n_4).h4)); 
L(n_4).Pvv_pJ_con   = ((-2)*(((L(n_4).L1vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvv_pJ       = (-2)*((-2)*(((L(n_4).L2vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvv_pJ_con   = (-2)*((-2)*(((L(n_4).L2vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Svv_pJ           = (-1i)*(sqrt(3)*k_m)*((((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Svv_pJ_con       = (-1i)*(sqrt(3)*k_p)*((((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                           + ((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h5)));
L(n_4).Svv_pJ_dag_con 	= (-1i)*(sqrt(3)*k_m)*((((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Svv_pJ_dag       = (-1i)*(sqrt(3)*k_p)*((((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                           + ((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Rvv_pJ       = 0; 
L(n_4).Rvv_pJ_dag   = 0;
L(n_4).Cvv_pJ       = 0;
L(n_4).Cvv_pJ_dag   = 0;
L(n_4).Zvv_pJ       = 0;        
L(n_4).Zvv_pJ_dag   = 0; 
L(n_4).Qvs_pJ       = (-2)*((-2)*(((L(n_4).L2vs)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvs_pJ_dag   = (-2)*((-2)*(((L(n_4).L2vs)*L(n_4).h)/L(n_4).h4));
L(n_4).Svs_pJ_dag_con  = (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h5))); 
L(n_4).Svs_pJ_dag       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g8m1_vs)/(L(n_4).h5))); 
L(n_4).Svs_pJ       	= (-1i)*(sqrt(6)*k_m)*(((L(n_4).Z_g8m1_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).Svs_pJ_con       = (-1i)*(sqrt(6)*k_p)*(((L(n_4).Z_g8m1_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_pJ_dag_con  = (-1i)*(sqrt(2/3)*k_m)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_pJ_dag     = (-1i)*(sqrt(2/3)*k_p)*(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_pJ         = (-1i)*(sqrt(2/3)*k_m)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_pJ_con     = (-1i)*(sqrt(2/3)*k_p)*(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                           + ((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).Rvs_pJ           = 0;
L(n_4).Rvs_pJ_dag_con   = 0;
L(n_4).Rvs_pJ_con       = 0;
L(n_4).Rvs_pJ_dag       = 0;
L(n_4).Ps_pJ        = ((-2)*((((L(n_4).L1ss))*L(n_4).h)/L(n_4).h4)); 
L(n_4).Ps_pJ_con    = ((-2)*((((L(n_4).L1ss))*L(n_4).h)/L(n_4).h4));
L(n_4).Cs_pJ        = 0;
L(n_4).Cs_pJ_dag    = 0;
%%%% Stage 3 Backward F{j=i-1} points
%%%% L(n_4).Aij_pI  =   +(-2)*((((gi-1) - ((gi-1) - (gi))/4)*L(n_4-1).h)/L(n_4).h1)    - ((gi) + (gi-1))/(2*L(n_4).h2) - ((gi-1)/(L(n_4).h2)); 
L(n_4).Pe_mJ        = ((-2)*((((L(n_4).L1eg7m))*L(n_4).h)/L(n_4).h4)); 
L(n_4).Pe_mJ_con    = ((-2)*((((L(n_4).L1eg7m))*L(n_4).h)/L(n_4).h4)); 
L(n_4).Ce_mJ        = 0;
L(n_4).Ce_mJ_dag    = 0;
L(n_4).Te_mJ       = 0;
L(n_4).Te_mJ_dag   = 0;
L(n_4).Ue_mJ       = (-1i)*(2)*-(1/2)*((((L(n_4).Xi_cv))/(L(n_4).h5)) + ((L(n_4).Xi_cv)/(L(n_4).h5))); 
L(n_4).We_mJ       = (-1i)*(1)*-(1/2)*((((L(n_4).Xi_cs))/(L(n_4).h5)) + ((L(n_4).Xi_cs)/(L(n_4).h5))); 
L(n_4).Ve_mJ       = 0;
L(n_4).Ve_mJ_dag   = 0;
L(n_4).Pvv_mJ       = ((-2)*(((L(n_4).L1vv)*L(n_4).h)/L(n_4).h4)); 
L(n_4).Pvv_mJ_con   = ((-2)*(((L(n_4).L1vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvv_mJ       = (-2)*((-2)*(((L(n_4).L2vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvv_mJ_con   = (-2)*((-2)*(((L(n_4).L2vv)*L(n_4).h)/L(n_4).h4));
L(n_4).Svv_mJ           = (-1i)*(sqrt(3)*k_m)*-((((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                             +((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Svv_mJ_con       = (-1i)*(sqrt(3)*k_p)*-((((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                             +((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv)/(L(n_4).h5)));
L(n_4).Svv_mJ_dag_con 	= (-1i)*(sqrt(3)*k_m)*-((((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                             +((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Svv_mJ_dag       = (-1i)*(sqrt(3)*k_p)*-((((2*L(n_4).Z_g7m_vv - 4*L(n_4).Z_g8m2_vv - L(n_4).Z_g8m3_vv))/(L(n_4).h5)) ...
                               +((2*L(n_4).Z_g6m_vv + L(n_4).Z_g8m3_vv)/(L(n_4).h5))); 
L(n_4).Rvv_mJ       = 0; 
L(n_4).Rvv_mJ_dag   = 0;
L(n_4).Cvv_mJ       = 0;
L(n_4).Cvv_mJ_dag   = 0;
L(n_4).Zvv_mJ       = 0;        
L(n_4).Zvv_mJ_dag   = 0; 
L(n_4).Qvs_mJ       = (-2)*((-2)*(((L(n_4).L2vs)*L(n_4).h)/L(n_4).h4));
L(n_4).Qvs_mJ_dag   = (-2)*((-2)*(((L(n_4).L2vs)*L(n_4).h)/L(n_4).h4));
L(n_4).Svs_mJ_dag_con   = (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                             +((L(n_4).Z_g8m1_vs)/(L(n_4).h5))); 
L(n_4).Svs_mJ_dag       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                             +((L(n_4).Z_g8m1_vs)/(L(n_4).h5))); 
L(n_4).Svs_mJ         	= (-1i)*(sqrt(6)*k_m)*-(((L(n_4).Z_g8m1_vs))/(L(n_4).h5) ...
                             +((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).Svs_mJ_con       = (-1i)*(sqrt(6)*k_p)*-(((L(n_4).Z_g8m1_vs))/(L(n_4).h5) ...
                               +((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_mJ_dag_con  = (-1i)*(sqrt(2/3)*k_m)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                               +((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_mJ_dag     = (-1i)*(sqrt(2/3)*k_p)*-(((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                               +((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_mJ     	= (-1i)*(sqrt(2/3)*k_m)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                               +((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).SIGvs_mJ_con     = (-1i)*(sqrt(2/3)*k_p)*-(((2*L(n_4).Z_g7m_vs + L(n_4).Z_g8m1_vs + 4*L(n_4).Z_g8m2_vs))/(L(n_4).h5) ...
                               +((L(n_4).Z_g7m_vs + 2*L(n_4).Z_g8m1_vs + 2*L(n_4).Z_g8m2_vs)/(L(n_4).h5))); 
L(n_4).Rvs_mJ       = 0;
L(n_4).Rvs_mJ_con   = 0;
L(n_4).Rvs_mJ_dag_con   = 0;
L(n_4).Rvs_mJ_dag   = 0;
L(n_4).Ps_mJ        = ((-2)*(((L(n_4).L1ss)*L(n_4).h)/L(n_4).h4)); 
L(n_4).Ps_mJ_con    = ((-2)*(((L(n_4).L1ss)*L(n_4).h)/L(n_4).h4));
L(n_4).Cs_mJ        = 0;
L(n_4).Cs_mJ_dag    = 0;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set up the Hamiltonian; H = H_E (Energy) + H_Epsilon (Strain) + 
%%%% K (First order k) + L (Second order k), for the conduction and valence
%%%% band. Each block is defined, "c" conduction band, "v" valence band "s"
%%%% spin orbit band.
%%%% Stage 1. (H_E (Energy)+ H_Epsilon (Strain)) blocks
L(n_4).E_centerI =  [(L(n_4).E_ccI + L(n_4).E_cI-Offset) (0)                           (0)  (0)  (0)  (0)  (0)  (0) 
                     (0)                          (L(n_4).E_ccI + L(n_4).E_cI-Offset)  (0)  (0)  (0)  (0)  (0)  (0)
                     (0)                          (0)                           (L(n_4).E_v1I + L(n_4).E_pI + L(n_4).E_qI-Offset) (-L(n_4).E_sI)                              (L(n_4).E_rI)                               (0)                                           (sqrt(3/2)*L(n_4).E_sI')        (sqrt(2)*L(n_4).E_qI) 
                     (0)                          (0)                          	(-L(n_4).E_sI')                            (L(n_4).E_v2I + L(n_4).E_pI - L(n_4).E_qI-Offset) 	(0)                                         (L(n_4).E_rI)                                 (sqrt(2)*L(n_4).E_rI')          (sqrt(1/2)*L(n_4).E_sI)
                     (0)                          (0)                         	(L(n_4).E_rI')                             (0)                                         (L(n_4).E_v2I + L(n_4).E_pI - L(n_4).E_qI-Offset)  (L(n_4).E_sI)                                 (sqrt(1/2)*L(n_4).E_sI)         (-sqrt(2)*L(n_4).E_rI)  
                     (0)                          (0)                           (0)                                        (L(n_4).E_rI')                              (L(n_4).E_sI')                              (L(n_4).E_v1I + L(n_4).E_pI + L(n_4).E_qI-Offset)    (-sqrt(2)*L(n_4).E_qI)          (sqrt(3/2)*L(n_4).E_sI)
                     (0)                          (0)                           (sqrt(3/2)*L(n_4).E_sI)                    (sqrt(2)*L(n_4).E_rI)                       (sqrt(1/2)*L(n_4).E_sI')                    (-sqrt(2)*L(n_4).E_qI)                        (L(n_4).E_ssI + L(n_4).E_pI-Offset)    (0)     
                     (0)                          (0)                           (sqrt(2)*L(n_4).E_qI')                     (sqrt(1/2)*L(n_4).E_sI)                     (-sqrt(2)*L(n_4).E_rI')                     (sqrt(3/2)*L(n_4).E_sI')                      (0)                             (L(n_4).E_ssI + L(n_4).E_pI-Offset)];
L(n_4).E_centerIR =  [(L(n_4).E_ccIR + L(n_4).E_cIR-Offset) (0)                           (0)  (0)  (0)  (0)  (0)  (0) 
                     (0)                          (L(n_4).E_ccIR + L(n_4).E_cIR-Offset)  (0)  (0)  (0)  (0)  (0)  (0)
                     (0)                          (0)                           (L(n_4).E_v1IR + L(n_4).E_pIR + L(n_4).E_qIR-Offset) (-L(n_4).E_sIR)                                 (L(n_4).E_rIR)                                  (0)                                             (sqrt(3/2)*L(n_4).E_sIR')   	(sqrt(2)*L(n_4).E_qIR) 
                     (0)                          (0)                           (-L(n_4).E_sIR')                               (L(n_4).E_v2IR + L(n_4).E_pIR - L(n_4).E_qIR-Offset) 	(0)                                             (L(n_4).E_rIR)                                  (sqrt(2)*L(n_4).E_rIR')     	(sqrt(1/2)*L(n_4).E_sIR)
                     (0)                          (0)                           (L(n_4).E_rIR')                                (0)                                             (L(n_4).E_v2IR + L(n_4).E_pIR - L(n_4).E_qIR-Offset)   (L(n_4).E_sIR)                                  (sqrt(1/2)*L(n_4).E_sIR)     	(-sqrt(2)*L(n_4).E_rIR)  
                     (0)                          (0)                           (0)                                            (L(n_4).E_rIR')                                 (L(n_4).E_sIR')                                 (L(n_4).E_v1IR + L(n_4).E_pIR + L(n_4).E_qIR-Offset) 	(-sqrt(2)*L(n_4).E_qIR)     	(sqrt(3/2)*L(n_4).E_sIR)
                     (0)                          (0)                           (sqrt(3/2)*L(n_4).E_sIR)                       (sqrt(2)*L(n_4).E_rIR)                          (sqrt(1/2)*L(n_4).E_sIR')                       (-sqrt(2)*L(n_4).E_qIR)                         (L(n_4).E_ssIR + L(n_4).E_pIR-Offset) 	(0)     
                     (0)                          (0)                           (sqrt(2)*L(n_4).E_qIR')                        (sqrt(1/2)*L(n_4).E_sIR)                        (-sqrt(2)*L(n_4).E_rIR')                        (sqrt(3/2)*L(n_4).E_sIR')                       (0)                             (L(n_4).E_ssIR + L(n_4).E_pIR-Offset)];
L(n_4).E_centerIL =  [(L(n_4).E_ccIL + L(n_4).E_cIL-Offset) (0)                           (0)  (0)  (0)  (0)  (0)  (0) 
                     (0)                          (L(n_4).E_ccIL + L(n_4).E_cIL-Offset)  (0)  (0)  (0)  (0)  (0)  (0)
                     (0)                          (0)                           (L(n_4).E_v1IL + L(n_4).E_pIL + L(n_4).E_qIL-Offset) (-L(n_4).E_sIL)                                 (L(n_4).E_rIL)                                  (0)                                             (sqrt(3/2)*L(n_4).E_sIL')    	(sqrt(2)*L(n_4).E_qIL) 
                     (0)                          (0)                           (-L(n_4).E_sIL')                               (L(n_4).E_v2IL + L(n_4).E_pIL - L(n_4).E_qIL-Offset) 	(0)                                             (L(n_4).E_rIL)                                  (sqrt(2)*L(n_4).E_rIL')      	(sqrt(1/2)*L(n_4).E_sIL)
                     (0)                          (0)                           (L(n_4).E_rIL')                                (0)                                             (L(n_4).E_v2IL + L(n_4).E_pIL - L(n_4).E_qIL-Offset)   (L(n_4).E_sIL)                                  (sqrt(1/2)*L(n_4).E_sIL)        (-sqrt(2)*L(n_4).E_rIL)  
                     (0)                          (0)                           (0)                                            (L(n_4).E_rIL')                                 (L(n_4).E_sIL')                                 (L(n_4).E_v1IL + L(n_4).E_pIL + L(n_4).E_qIL-Offset)  	(-sqrt(2)*L(n_4).E_qIL)      	(sqrt(3/2)*L(n_4).E_sIL)
                     (0)                          (0)                           (sqrt(3/2)*L(n_4).E_sIL)                       (sqrt(2)*L(n_4).E_rIL)                          (sqrt(1/2)*L(n_4).E_sIL')                       (-sqrt(2)*L(n_4).E_qIL)                         (L(n_4).E_ssIL + L(n_4).E_pIL-Offset)	(0)     
                     (0)                          (0)                           (sqrt(2)*L(n_4).E_qIL')                        (sqrt(1/2)*L(n_4).E_sIL)                        (-sqrt(2)*L(n_4).E_rIL')                        (sqrt(3/2)*L(n_4).E_sIL')                       (0)                             (L(n_4).E_ssIL + L(n_4).E_pIL-Offset)];
L(n_4).E_centerJ =  [(L(n_4).E_ccJ + L(n_4).E_cJ-Offset) (0)                           (0)  (0)  (0)  (0)  (0)  (0) 
                     (0)                          (L(n_4).E_ccJ + L(n_4).E_cJ-Offset)  (0)  (0)  (0)  (0)  (0)  (0)
                     (0)                          (0)                           (L(n_4).E_v1J + L(n_4).E_pJ + L(n_4).E_qJ-Offset) (-L(n_4).E_sJ)                              (L(n_4).E_rJ)                               (0)                                           (sqrt(3/2)*L(n_4).E_sJ')        (sqrt(2)*L(n_4).E_qJ) 
                     (0)                          (0)                           (-L(n_4).E_sJ')                            (L(n_4).E_v2J + L(n_4).E_pJ - L(n_4).E_qJ-Offset) 	(0)                                         (L(n_4).E_rJ)                                 (sqrt(2)*L(n_4).E_rJ')          (sqrt(1/2)*L(n_4).E_sJ)
                     (0)                          (0)                           (L(n_4).E_rJ')                             (0)                                         (L(n_4).E_v2J + L(n_4).E_pJ - L(n_4).E_qJ-Offset)  (L(n_4).E_sJ)                                 (sqrt(1/2)*L(n_4).E_sJ)         (-sqrt(2)*L(n_4).E_rJ)  
                     (0)                          (0)                           (0)                                        (L(n_4).E_rJ')                              (L(n_4).E_sJ')                              (L(n_4).E_v1J + L(n_4).E_pJ + L(n_4).E_qJ-Offset)    (-sqrt(2)*L(n_4).E_qJ)          (sqrt(3/2)*L(n_4).E_sJ)
                     (0)                          (0)                           (sqrt(3/2)*L(n_4).E_sJ)                    (sqrt(2)*L(n_4).E_rJ)                       (sqrt(1/2)*L(n_4).E_sJ')                    (-sqrt(2)*L(n_4).E_qJ)                        (L(n_4).E_ssJ + L(n_4).E_pJ-Offset)    (0)     
                     (0)                          (0)                           (sqrt(2)*L(n_4).E_qJ')                     (sqrt(1/2)*L(n_4).E_sJ)                     (-sqrt(2)*L(n_4).E_rJ')                     (sqrt(3/2)*L(n_4).E_sJ')                      (0)                             (L(n_4).E_ssJ + L(n_4).E_pJ-Offset)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Stage 2. K (First order k) blocks
%%%% Interface Kcv
L(n_4).Kcv_centerI =  [(L(n_4).Te_cI) (0)                     (sqrt(3)*L(n_4).Te_cI_dag) (L(n_4).Ue_cI)
                       (L(n_4).Ue_cI) (-sqrt(3)*L(n_4).Te_cI) (0)                        (-L(n_4).Te_cI_dag) ];                  
L(n_4).Kcv_centerIR =  [(L(n_4).Te_cIR) (0)                     (sqrt(3)*L(n_4).Te_cIR_dag) (L(n_4).Ue_cIR)
                       (L(n_4).Ue_cIR) (-sqrt(3)*L(n_4).Te_cIR) (0)                        (-L(n_4).Te_cIR_dag) ];                  
L(n_4).Kcv_centerIL =  [(L(n_4).Te_cIL) (0)                     (sqrt(3)*L(n_4).Te_cIL_dag) (L(n_4).Ue_cIL)
                       (L(n_4).Ue_cIL) (-sqrt(3)*L(n_4).Te_cIL) (0)                        (-L(n_4).Te_cIL_dag) ];
L(n_4).Kcv_plusI =  [(L(n_4).Te_pI) (0)                     (sqrt(3)*L(n_4).Te_pI_dag) (L(n_4).Ue_pI)
                       (L(n_4).Ue_pI) (-sqrt(3)*L(n_4).Te_pI) (0)                        (-L(n_4).Te_pI_dag) ];
L(n_4).Kcv_minusI =  [(L(n_4).Te_mI) (0)                     (sqrt(3)*L(n_4).Te_mI_dag) (L(n_4).Ue_mI)
                       (L(n_4).Ue_mI) (-sqrt(3)*L(n_4).Te_mI) (0)                        (-L(n_4).Te_mI_dag) ];                  
%%%% Layer Kcv                   
L(n_4).Kcv_centerJ =  [(L(n_4).Te_cJ) (0)                     (sqrt(3)*L(n_4).Te_cJ_dag) (L(n_4).Ue_cJ)
                       (L(n_4).Ue_cJ) (-sqrt(3)*L(n_4).Te_cJ) (0)                        (-L(n_4).Te_cJ_dag) ];
L(n_4).Kcv_plusJ =  [(L(n_4).Te_pJ) (0)                     (sqrt(3)*L(n_4).Te_pJ_dag) (L(n_4).Ue_pJ)
                       (L(n_4).Ue_pJ) (-sqrt(3)*L(n_4).Te_pJ) (0)                        (-L(n_4).Te_pJ_dag) ];
L(n_4).Kcv_minusJ =  [(L(n_4).Te_mJ) (0)                     (sqrt(3)*L(n_4).Te_mJ_dag) (L(n_4).Ue_mJ)
                       (L(n_4).Ue_mJ) (-sqrt(3)*L(n_4).Te_mJ) (0)                        (-L(n_4).Te_mJ_dag) ];                                   
%%%% Interface Kvc
L(n_4).Kvc_centerI =  [(L(n_4).Te_cI_dag)   	(L(n_4).Ue_cI)
                       (0)                      (-sqrt(3)*L(n_4).Te_cI_dag)
                       (sqrt(3)*L(n_4).Te_cI)   (0) 
                       (L(n_4).Ue_cI)         	(-L(n_4).Te_cI) ];
L(n_4).Kvc_centerIR =  [(L(n_4).Te_cIR_dag)   	(L(n_4).Ue_cIR)
                       (0)                      (-sqrt(3)*L(n_4).Te_cIR_dag)
                       (sqrt(3)*L(n_4).Te_cIR)   (0) 
                       (L(n_4).Ue_cIR)         	(-L(n_4).Te_cIR) ];
L(n_4).Kvc_centerIL =  [(L(n_4).Te_cIL_dag)   	(L(n_4).Ue_cIL)
                       (0)                      (-sqrt(3)*L(n_4).Te_cIL_dag)
                       (sqrt(3)*L(n_4).Te_cIL)   (0) 
                       (L(n_4).Ue_cIL)         	(-L(n_4).Te_cIL) ];                                     
L(n_4).Kvc_plusI    =  [(L(n_4).Te_pI_dag)   	(L(n_4).Ue_pI)
                       (0)                      (-sqrt(3)*L(n_4).Te_pI_dag)
                       (sqrt(3)*L(n_4).Te_pI)   (0) 
                       (L(n_4).Ue_pI)         	(-L(n_4).Te_pI) ];                    
L(n_4).Kvc_minusI  =  [(L(n_4).Te_mI_dag)   	(L(n_4).Ue_mI)
                       (0)                      (-sqrt(3)*L(n_4).Te_mI_dag)
                       (sqrt(3)*L(n_4).Te_mI)   (0) 
                       (L(n_4).Ue_mI)         	(-L(n_4).Te_mI) ];                    
%%%% Layer Kvc                    
L(n_4).Kvc_centerJ =  [(L(n_4).Te_cJ_dag)   	(L(n_4).Ue_cJ)
                       (0)                      (-sqrt(3)*L(n_4).Te_cJ_dag)
                       (sqrt(3)*L(n_4).Te_cJ)   (0) 
                       (L(n_4).Ue_cJ)         	(-L(n_4).Te_cJ) ];
L(n_4).Kvc_plusJ    =  [(L(n_4).Te_pJ_dag)   	(L(n_4).Ue_pJ)
                       (0)                      (-sqrt(3)*L(n_4).Te_pJ_dag)
                       (sqrt(3)*L(n_4).Te_pJ)   (0) 
                       (L(n_4).Ue_pJ)         	(-L(n_4).Te_pJ) ];                    
L(n_4).Kvc_minusJ  =  [(L(n_4).Te_mJ_dag)   	(L(n_4).Ue_mJ)
                       (0)                      (-sqrt(3)*L(n_4).Te_mJ_dag)
                       (sqrt(3)*L(n_4).Te_mJ)   (0) 
                       (L(n_4).Ue_mJ)         	(-L(n_4).Te_mJ) ]; 
%%%% Interface Kcs                   
L(n_4).Kcs_centerI =  [(L(n_4).We_cI)     (L(n_4).Ve_cI) 
                       (L(n_4).Ve_cI_dag) (-L(n_4).We_cI) ];                   
L(n_4).Kcs_centerIR =  [(L(n_4).We_cIR)     (L(n_4).Ve_cIR) 
                       (L(n_4).Ve_cIR_dag) (-L(n_4).We_cIR) ];                     
L(n_4).Kcs_centerIL =  [(L(n_4).We_cIL)     (L(n_4).Ve_cIL) 
                       (L(n_4).Ve_cIL_dag) (-L(n_4).We_cIL) ];                    
L(n_4).Kcs_plusI =  [(L(n_4).We_pI)     (L(n_4).Ve_pI) 
                       (L(n_4).Ve_pI_dag) (-L(n_4).We_pI) ];                  
L(n_4).Kcs_minusI =  [(L(n_4).We_mI)     (L(n_4).Ve_mI) 
                       (L(n_4).Ve_mI_dag) (-L(n_4).We_mI) ];                   
%%%% Layer Kcs                      
L(n_4).Kcs_centerJ =  [(L(n_4).We_cJ)     (L(n_4).Ve_cJ) 
                       (L(n_4).Ve_cJ_dag) (-L(n_4).We_cJ) ];                                      
L(n_4).Kcs_plusJ =  [(L(n_4).We_pJ)     (L(n_4).Ve_pJ) 
                       (L(n_4).Ve_pJ_dag) (-L(n_4).We_pJ) ];                   
L(n_4).Kcs_minusJ =  [(L(n_4).We_mJ)     (L(n_4).Ve_mJ) 
                       (L(n_4).Ve_mJ_dag) (-L(n_4).We_mJ) ];                                       
%%%% Interface Ksc                   
L(n_4).Ksc_centerI =  [(L(n_4).We_cI)     (L(n_4).Ve_cI) 
                       (L(n_4).Ve_cI_dag) (-L(n_4).We_cI) ];                  
L(n_4).Ksc_centerIR =  [(L(n_4).We_cIR)     (L(n_4).Ve_cIR) 
                       (L(n_4).Ve_cIR_dag) (-L(n_4).We_cIR) ];                     
L(n_4).Ksc_centerIL =  [(L(n_4).We_cIL)     (L(n_4).Ve_cIL) 
                       (L(n_4).Ve_cIL_dag) (-L(n_4).We_cIL) ];                    
L(n_4).Ksc_plusI =  [(L(n_4).We_pI)     (L(n_4).Ve_pI) 
                       (L(n_4).Ve_pI_dag) (-L(n_4).We_pI) ];                   
L(n_4).Ksc_minusI =  [(L(n_4).We_mI)     (L(n_4).Ve_mI) 
                       (L(n_4).Ve_mI_dag) (-L(n_4).We_mI) ];                   
%%%% Layer Ksc                      
L(n_4).Ksc_centerJ =  [(L(n_4).We_cJ)     (L(n_4).Ve_cJ) 
                       (L(n_4).Ve_cJ_dag) (-L(n_4).We_cJ) ];                                     
L(n_4).Ksc_plusJ =  [(L(n_4).We_pJ)     (L(n_4).Ve_pJ) 
                       (L(n_4).Ve_pJ_dag) (-L(n_4).We_pJ) ];                   
L(n_4).Ksc_minusJ =  [(L(n_4).We_mJ)     (L(n_4).Ve_mJ) 
                       (L(n_4).Ve_mJ_dag) (-L(n_4).We_mJ) ];                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Stage 3. L (Second order k) blocks
%%%% Interface Lcc 
L(n_4).Lcc_centerI =  [(L(n_4).Pe_cI)     (L(n_4).Ce_cI) 
                       (L(n_4).Ce_cI_dag) (L(n_4).Pe_cI_con) ];
L(n_4).Lcc_centerIR =  [(L(n_4).Pe_cIR)     (L(n_4).Ce_cIR) 
                       (L(n_4).Ce_cIR_dag)  (L(n_4).Pe_cIR_con) ];                                      
L(n_4).Lcc_centerIL =  [(L(n_4).Pe_cIL)     (L(n_4).Ce_cIL) 
                       (L(n_4).Ce_cIL_dag)  (L(n_4).Pe_cIL_con) ];                   
L(n_4).Lcc_plusI   =  [(L(n_4).Pe_pI)     (L(n_4).Ce_pI) 
                       (L(n_4).Ce_pI_dag) (L(n_4).Pe_pI_con) ];
L(n_4).Lcc_minusI  =  [(L(n_4).Pe_mI)     (L(n_4).Ce_mI) 
                       (L(n_4).Ce_mI_dag) (L(n_4).Pe_mI_con) ];
%%%% Layer Lcc                    
L(n_4).Lcc_centerJ =  [(L(n_4).Pe_cJ)     (L(n_4).Ce_cJ) 
                       (L(n_4).Ce_cJ_dag) (L(n_4).Pe_cJ_con) ];
L(n_4).Lcc_plusJ   =  [(L(n_4).Pe_pJ)     (L(n_4).Ce_pJ) 
                       (L(n_4).Ce_pJ_dag) (L(n_4).Pe_pJ_con) ];
L(n_4).Lcc_minusJ  =  [(L(n_4).Pe_mJ)     (L(n_4).Ce_mJ) 
                       (L(n_4).Ce_mJ_dag) (L(n_4).Pe_mJ_con) ];   
%%%% Interface Lvv 
L(n_4).Lvv_centerI =  [(L(n_4).Pvv_cI + L(n_4).Qvv_cI) (-L(n_4).Svv_cI)                (L(n_4).Rvv_cI)                         (L(n_4).Cvv_cI_dag)
                       (-L(n_4).Svv_cI_dag)            (L(n_4).Pvv_cI - L(n_4).Qvv_cI) (L(n_4).Zvv_cI)                         (L(n_4).Rvv_cI)
                       (L(n_4).Rvv_cI_dag)             (L(n_4).Zvv_cI_dag)             (L(n_4).Pvv_cI_con - L(n_4).Qvv_cI_con) (L(n_4).Svv_cI_dag_con)
                       (L(n_4).Cvv_cI)                 (L(n_4).Rvv_cI_dag)             (L(n_4).Svv_cI_con)                     (L(n_4).Pvv_cI_con + L(n_4).Qvv_cI_con) ];
L(n_4).Lvv_centerIR =  [(L(n_4).Pvv_cIR + L(n_4).Qvv_cIR) (-L(n_4).Svv_cIR)                 (L(n_4).Rvv_cIR)                          (L(n_4).Cvv_cIR_dag)
                       (-L(n_4).Svv_cIR_dag)              (L(n_4).Pvv_cIR - L(n_4).Qvv_cIR) (L(n_4).Zvv_cIR)                          (L(n_4).Rvv_cIR)
                       (L(n_4).Rvv_cIR_dag)               (L(n_4).Zvv_cIR_dag)              (L(n_4).Pvv_cIR_con - L(n_4).Qvv_cIR_con) (L(n_4).Svv_cIR_dag_con)
                       (L(n_4).Cvv_cIR)                   (L(n_4).Rvv_cIR_dag)              (L(n_4).Svv_cIR_con)                      (L(n_4).Pvv_cIR_con + L(n_4).Qvv_cIR_con) ];                   
L(n_4).Lvv_centerIL =  [(L(n_4).Pvv_cIL + L(n_4).Qvv_cIL) (-L(n_4).Svv_cIL)                 (L(n_4).Rvv_cIL)                          (L(n_4).Cvv_cIL_dag)
                       (-L(n_4).Svv_cIL_dag)              (L(n_4).Pvv_cIL - L(n_4).Qvv_cIL) (L(n_4).Zvv_cIL)                          (L(n_4).Rvv_cIL)
                       (L(n_4).Rvv_cIL_dag)               (L(n_4).Zvv_cIL_dag)              (L(n_4).Pvv_cIL_con - L(n_4).Qvv_cIL_con) (L(n_4).Svv_cIL_dag_con)
                       (L(n_4).Cvv_cIL)                   (L(n_4).Rvv_cIL_dag)              (L(n_4).Svv_cIL_con)                      (L(n_4).Pvv_cIL_con + L(n_4).Qvv_cIL_con) ];                                     
L(n_4).Lvv_plusI =    [(L(n_4).Pvv_pI + L(n_4).Qvv_pI) (-L(n_4).Svv_pI)                (L(n_4).Rvv_pI)                         (L(n_4).Cvv_pI_dag)
                       (-L(n_4).Svv_pI_dag)            (L(n_4).Pvv_pI - L(n_4).Qvv_pI) (L(n_4).Zvv_pI)                         (L(n_4).Rvv_pI)
                       (L(n_4).Rvv_pI_dag)             (L(n_4).Zvv_pI_dag)             (L(n_4).Pvv_pI_con - L(n_4).Qvv_pI_con) (L(n_4).Svv_pI_dag_con)
                       (L(n_4).Cvv_pI)                 (L(n_4).Rvv_pI_dag)             (L(n_4).Svv_pI_con)                     (L(n_4).Pvv_pI_con + L(n_4).Qvv_pI_con) ];
L(n_4).Lvv_minusI =   [(L(n_4).Pvv_mI + L(n_4).Qvv_mI) (-L(n_4).Svv_mI)                (L(n_4).Rvv_mI)                         (L(n_4).Cvv_mI_dag)
                       (-L(n_4).Svv_mI_dag)            (L(n_4).Pvv_mI - L(n_4).Qvv_mI) (L(n_4).Zvv_mI)                         (L(n_4).Rvv_mI)
                       (L(n_4).Rvv_mI_dag)             (L(n_4).Zvv_mI_dag)             (L(n_4).Pvv_mI_con - L(n_4).Qvv_mI_con) (L(n_4).Svv_mI_dag_con)
                       (L(n_4).Cvv_mI)                 (L(n_4).Rvv_mI_dag)             (L(n_4).Svv_mI_con)                     (L(n_4).Pvv_mI_con + L(n_4).Qvv_mI_con) ];                   
%%%% Layer Lvv                    
L(n_4).Lvv_centerJ =  [(L(n_4).Pvv_cJ + L(n_4).Qvv_cJ) (-L(n_4).Svv_cJ)                (L(n_4).Rvv_cJ)                         (L(n_4).Cvv_cJ_dag)
                       (-L(n_4).Svv_cJ_dag)            (L(n_4).Pvv_cJ - L(n_4).Qvv_cJ) (L(n_4).Zvv_cJ)                         (L(n_4).Rvv_cJ)
                       (L(n_4).Rvv_cJ_dag)             (L(n_4).Zvv_cJ_dag)             (L(n_4).Pvv_cJ_con - L(n_4).Qvv_cJ_con) (L(n_4).Svv_cJ_dag_con)
                       (L(n_4).Cvv_cJ)                 (L(n_4).Rvv_cJ_dag)             (L(n_4).Svv_cJ_con)                     (L(n_4).Pvv_cJ_con + L(n_4).Qvv_cJ_con) ];
L(n_4).Lvv_plusJ =    [(L(n_4).Pvv_pJ + L(n_4).Qvv_pJ) (-L(n_4).Svv_pJ)                (L(n_4).Rvv_pJ)                         (L(n_4).Cvv_pJ_dag)
                       (-L(n_4).Svv_pJ_dag)            (L(n_4).Pvv_pJ - L(n_4).Qvv_pJ) (L(n_4).Zvv_pJ)                         (L(n_4).Rvv_pJ)
                       (L(n_4).Rvv_pJ_dag)             (L(n_4).Zvv_pJ_dag)             (L(n_4).Pvv_pJ_con - L(n_4).Qvv_pJ_con) (L(n_4).Svv_pJ_dag_con)
                       (L(n_4).Cvv_pJ)                 (L(n_4).Rvv_pJ_dag)             (L(n_4).Svv_pJ_con)                     (L(n_4).Pvv_pJ_con + L(n_4).Qvv_pJ_con) ];
L(n_4).Lvv_minusJ =   [(L(n_4).Pvv_mJ + L(n_4).Qvv_mJ) (-L(n_4).Svv_mJ)                (L(n_4).Rvv_mJ)                         (L(n_4).Cvv_mJ_dag)
                       (-L(n_4).Svv_mJ_dag)            (L(n_4).Pvv_mJ - L(n_4).Qvv_mJ) (L(n_4).Zvv_mJ)                         (L(n_4).Rvv_mJ)
                       (L(n_4).Rvv_mJ_dag)             (L(n_4).Zvv_mJ_dag)             (L(n_4).Pvv_mJ_con - L(n_4).Qvv_mJ_con) (L(n_4).Svv_mJ_dag_con)
                       (L(n_4).Cvv_mJ)                 (L(n_4).Rvv_mJ_dag)             (L(n_4).Svv_mJ_con)                     (L(n_4).Pvv_mJ_con + L(n_4).Qvv_mJ_con) ];       
%%%% Interface Lvs
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
L(n_4).Lvs_minusI = [(sqrt(3/2)*L(n_4).SIGvs_mI_con) (sqrt(2)*L(n_4).Qvs_mI)                 
                       (sqrt(2)*L(n_4).Rvs_mI_con)     (sqrt(1/2)*L(n_4).Svs_mI_con)
                       (sqrt(1/2)*L(n_4).Svs_mI)       (-sqrt(2)*L(n_4).Rvs_mI)         
                       (-sqrt(2)*L(n_4).Qvs_mI)        (sqrt(3/2)*L(n_4).SIGvs_mI)  ];                   
%%%% Layer Lvs                     
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
%%%% Interface Lsv
L(n_4).Lsv_centerI  =  [(sqrt(3/2)*L(n_4).SIGvs_cI_dag_con) (sqrt(2)*L(n_4).Rvs_cI_dag_con)   (sqrt(1/2)*L(n_4).Svs_cI_dag) (-sqrt(2)*L(n_4).Qvs_cI_dag)         
                        (sqrt(2)*L(n_4).Qvs_cI_dag)         (sqrt(1/2)*L(n_4).Svs_cI_dag_con) (-sqrt(2)*L(n_4).Rvs_cI_dag)  (sqrt(3/2)*L(n_4).SIGvs_cI_dag)];                   
L(n_4).Lsv_centerIR =  [(sqrt(3/2)*L(n_4).SIGvs_cIR_dag_con) (sqrt(2)*L(n_4).Rvs_cIR_dag_con)   (sqrt(1/2)*L(n_4).Svs_cIR_dag) (-sqrt(2)*L(n_4).Qvs_cIR_dag)         
                        (sqrt(2)*L(n_4).Qvs_cIR_dag)         (sqrt(1/2)*L(n_4).Svs_cIR_dag_con) (-sqrt(2)*L(n_4).Rvs_cIR_dag)  (sqrt(3/2)*L(n_4).SIGvs_cIR_dag)];                   
L(n_4).Lsv_centerIL  =  [(sqrt(3/2)*L(n_4).SIGvs_cIL_dag_con) (sqrt(2)*L(n_4).Rvs_cIL_dag_con)   (sqrt(1/2)*L(n_4).Svs_cIL_dag) (-sqrt(2)*L(n_4).Qvs_cIL_dag)         
                        (sqrt(2)*L(n_4).Qvs_cIL_dag)         (sqrt(1/2)*L(n_4).Svs_cIL_dag_con) (-sqrt(2)*L(n_4).Rvs_cIL_dag)  (sqrt(3/2)*L(n_4).SIGvs_cIL_dag)];
L(n_4).Lsv_plusI    =  [(sqrt(3/2)*L(n_4).SIGvs_pI_dag_con) (sqrt(2)*L(n_4).Rvs_pI_dag_con)   (sqrt(1/2)*L(n_4).Svs_pI_dag) (-sqrt(2)*L(n_4).Qvs_pI_dag)         
                        (sqrt(2)*L(n_4).Qvs_pI_dag)         (sqrt(1/2)*L(n_4).Svs_pI_dag_con) (-sqrt(2)*L(n_4).Rvs_pI_dag)  (sqrt(3/2)*L(n_4).SIGvs_pI_dag)];
L(n_4).Lsv_minusI   =  [(sqrt(3/2)*L(n_4).SIGvs_mI_dag_con) (sqrt(2)*L(n_4).Rvs_mI_dag_con)   (sqrt(1/2)*L(n_4).Svs_mI_dag) (-sqrt(2)*L(n_4).Qvs_mI_dag)         
                        (sqrt(2)*L(n_4).Qvs_mI_dag)         (sqrt(1/2)*L(n_4).Svs_mI_dag_con) (-sqrt(2)*L(n_4).Rvs_mI_dag)  (sqrt(3/2)*L(n_4).SIGvs_mI_dag)];                   
%%%% Layer Lvs                    
L(n_4).Lsv_centerJ  =  [(sqrt(3/2)*L(n_4).SIGvs_cJ_dag_con) (sqrt(2)*L(n_4).Rvs_cJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_cJ_dag) (-sqrt(2)*L(n_4).Qvs_cJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_cJ_dag)         (sqrt(1/2)*L(n_4).Svs_cJ_dag_con) (-sqrt(2)*L(n_4).Rvs_cJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_cJ_dag)];
L(n_4).Lsv_plusJ    =  [(sqrt(3/2)*L(n_4).SIGvs_pJ_dag_con) (sqrt(2)*L(n_4).Rvs_pJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_pJ_dag) (-sqrt(2)*L(n_4).Qvs_pJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_pJ_dag)         (sqrt(1/2)*L(n_4).Svs_pJ_dag_con) (-sqrt(2)*L(n_4).Rvs_pJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_pJ_dag)];
L(n_4).Lsv_minusJ   =  [(sqrt(3/2)*L(n_4).SIGvs_mJ_dag_con) (sqrt(2)*L(n_4).Rvs_mJ_dag_con)   (sqrt(1/2)*L(n_4).Svs_mJ_dag) (-sqrt(2)*L(n_4).Qvs_mJ_dag)         
                        (sqrt(2)*L(n_4).Qvs_mJ_dag)         (sqrt(1/2)*L(n_4).Svs_mJ_dag_con) (-sqrt(2)*L(n_4).Rvs_mJ_dag)  (sqrt(3/2)*L(n_4).SIGvs_mJ_dag)];
%%%% Interface Lss
L(n_4).Lss_centerI  =  [(L(n_4).Ps_cI)     (L(n_4).Cs_cI) 
                       (L(n_4).Cs_cI_dag) (L(n_4).Ps_cI_con) ];
L(n_4).Lss_centerIR =  [(L(n_4).Ps_cIR)     (L(n_4).Cs_cIR) 
                       (L(n_4).Cs_cIR_dag)  (L(n_4).Ps_cIR_con) ];                                      
L(n_4).Lss_centerIL =  [(L(n_4).Ps_cIL)     (L(n_4).Cs_cIL) 
                       (L(n_4).Cs_cIL_dag)  (L(n_4).Ps_cIL_con) ];                   
L(n_4).Lss_plusI    =  [(L(n_4).Ps_pI)     (L(n_4).Cs_pI) 
                       (L(n_4).Cs_pI_dag) (L(n_4).Ps_pI_con) ];
L(n_4).Lss_minusI   =  [(L(n_4).Ps_mI)     (L(n_4).Cs_mI) 
                       (L(n_4).Cs_mI_dag) (L(n_4).Ps_mI_con) ];
%%%% Layer Lss                   
L(n_4).Lss_centerJ =  [(L(n_4).Ps_cJ)     (L(n_4).Cs_cJ) 
                       (L(n_4).Cs_cJ_dag) (L(n_4).Ps_cJ_con) ];
L(n_4).Lss_plusJ   =  [(L(n_4).Ps_pJ)     (L(n_4).Cs_pJ) 
                       (L(n_4).Cs_pJ_dag) (L(n_4).Ps_pJ_con) ];
L(n_4).Lss_minusJ  =  [(L(n_4).Ps_mJ)     (L(n_4).Cs_mJ) 
                       (L(n_4).Cs_mJ_dag) (L(n_4).Ps_mJ_con) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Stage 4. Combine each interaction to define Hamiltonian
L(n_4).Hkp_centerI = [ (zeros(2,2))          (L(n_4).Kcv_centerI)  (L(n_4).Kcs_centerI)
                       (L(n_4).Kvc_centerI)	 (zeros(4,4))          (zeros(4,2))
                       (L(n_4).Ksc_centerI)	 (zeros(2,4))          (zeros(2,2)) ];
L(n_4).Hkp_centerIR= [ (zeros(2,2))          (L(n_4).Kcv_centerIR) (L(n_4).Kcs_centerIR)
                       (L(n_4).Kvc_centerIR) (zeros(4,4))          (zeros(4,2))
                       (L(n_4).Ksc_centerIR) (zeros(2,4))          (zeros(2,2)) ];
L(n_4).Hkp_centerIL= [ (zeros(2,2))          (L(n_4).Kcv_centerIL) (L(n_4).Kcs_centerIL)
                       (L(n_4).Kvc_centerIL) (zeros(4,4))          (zeros(4,2))
                       (L(n_4).Ksc_centerIL) (zeros(2,4))          (zeros(2,2)) ];
L(n_4).Hkp_plusI   = [ (zeros(2,2))         (L(n_4).Kcv_plusI)	(L(n_4).Kcs_plusI)
                       (L(n_4).Kvc_plusI)	(zeros(4,4))        (zeros(4,2))
                       (L(n_4).Ksc_plusI)	(zeros(2,4))        (zeros(2,2)) ];                 
L(n_4).Hkp_minusI  = [ (zeros(2,2))          (L(n_4).Kcv_minusI) (L(n_4).Kcs_minusI)
                       (L(n_4).Kvc_minusI)	(zeros(4,4))         (zeros(4,2))
                       (L(n_4).Ksc_minusI)	(zeros(2,4))         (zeros(2,2)) ];                 
L(n_4).HL_centerI = [  (L(n_4).Lcc_centerI) (zeros(2,4))           (zeros(2,2))
                       (zeros(4,2))         ((L(n_4).Lvv_centerI)) (L(n_4).Lvs_centerI)
                       (zeros(2,2))         (L(n_4).Lsv_centerI)   ((L(n_4).Lss_centerI)) ];
L(n_4).HL_centerIR= [  (L(n_4).Lcc_centerIR) (zeros(2,4))            (zeros(2,2))
                       (zeros(4,2))          ((L(n_4).Lvv_centerIR)) (L(n_4).Lvs_centerIR)
                       (zeros(2,2))          (L(n_4).Lsv_centerIR)   ((L(n_4).Lss_centerIR)) ];
L(n_4).HL_centerIL= [  (L(n_4).Lcc_centerIL) (zeros(2,4))            (zeros(2,2))
                       (zeros(4,2))          ((L(n_4).Lvv_centerIL)) (L(n_4).Lvs_centerIL)
                       (zeros(2,2))          (L(n_4).Lsv_centerIL)   ((L(n_4).Lss_centerIL)) ];
L(n_4).HL_plusI   = [ (L(n_4).Lcc_plusI)   (zeros(2,4))                (zeros(2,2))
                      (zeros(4,2))         (L(n_4).Lvv_plusI) (L(n_4).Lvs_plusI)
                      (zeros(2,2))         (L(n_4).Lsv_plusI) (L(n_4).Lss_plusI) ];                 
L(n_4).HL_minusI  = [ (L(n_4).Lcc_minusI)  (zeros(2,4))        (zeros(2,2))
                      (zeros(4,2))         (L(n_4).Lvv_minusI) (L(n_4).Lvs_minusI)
                      (zeros(2,2))         (L(n_4).Lsv_minusI) (L(n_4).Lss_minusI) ];                
L(n_4).H_centerI   = (sqrt(Const).*L(n_4).Hkp_centerI  + -Const.*L(n_4).HL_centerI)  + L(n_4).E_centerI;
L(n_4).H_centerIR  = (sqrt(Const).*L(n_4).Hkp_centerIR + -Const.*L(n_4).HL_centerIR) + L(n_4).E_centerIR;
L(n_4).H_centerIL  = (sqrt(Const).*L(n_4).Hkp_centerIL + -Const.*L(n_4).HL_centerIL) + L(n_4).E_centerIL;
L(n_4).H_plusI     = (sqrt(Const).*L(n_4).Hkp_plusI    + -Const.*L(n_4).HL_plusI);
L(n_4).H_minusI    = (sqrt(Const).*L(n_4).Hkp_minusI   + -Const.*L(n_4).HL_minusI);                 
L(n_4).Hkp_centerJ = [ (zeros(2,2))         (L(n_4).Kcv_centerJ)	(L(n_4).Kcs_centerJ)
                       (L(n_4).Kvc_centerJ)	(zeros(4,4))            (zeros(4,2))
                       (L(n_4).Ksc_centerJ)	(zeros(2,4))            (zeros(2,2)) ];
L(n_4).Hkp_plusJ   = [ (zeros(2,2))         (L(n_4).Kcv_plusJ)	(L(n_4).Kcs_plusJ)
                       (L(n_4).Kvc_plusJ)	(zeros(4,4))        (zeros(4,2))
                       (L(n_4).Ksc_plusJ)	(zeros(2,4))        (zeros(2,2)) ];                
L(n_4).Hkp_minusJ  = [ (zeros(2,2))         (L(n_4).Kcv_minusJ)	(L(n_4).Kcs_minusJ)
                       (L(n_4).Kvc_minusJ)	(zeros(4,4))        (zeros(4,2))
                       (L(n_4).Ksc_minusJ)	(zeros(2,4))        (zeros(2,2)) ];                
L(n_4).HL_centerJ = [  (L(n_4).Lcc_centerJ) (zeros(2,4))           (zeros(2,2))
                       (zeros(4,2))         ((L(n_4).Lvv_centerJ)) (L(n_4).Lvs_centerJ)
                       (zeros(2,2))         (L(n_4).Lsv_centerJ)   ((L(n_4).Lss_centerJ)) ];
L(n_4).HL_plusJ   = [ (L(n_4).Lcc_plusJ)   (zeros(2,4))       (zeros(2,2))
                      (zeros(4,2))         (L(n_4).Lvv_plusJ) (L(n_4).Lvs_plusJ)
                      (zeros(2,2))         (L(n_4).Lsv_plusJ) (L(n_4).Lss_plusJ) ];                 
L(n_4).HL_minusJ  = [ (L(n_4).Lcc_minusJ)  (zeros(2,4))        (zeros(2,2))
                      (zeros(4,2))         (L(n_4).Lvv_minusJ) (L(n_4).Lvs_minusJ)
                      (zeros(2,2))         (L(n_4).Lsv_minusJ) (L(n_4).Lss_minusJ) ];                  
L(n_4).H_centerJ   = (sqrt(Const).*L(n_4).Hkp_centerJ + -Const.*L(n_4).HL_centerJ) + L(n_4).E_centerJ;
L(n_4).H_plusJ     = ((sqrt(Const).*L(n_4).Hkp_plusJ + -Const.*L(n_4).HL_plusJ));
L(n_4).H_minusJ    = ((sqrt(Const).*L(n_4).Hkp_minusJ + -Const.*L(n_4).HL_minusJ));
end      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Stage 2.4 Set up the block tridiagonal matrix by combining the 
%%%% digonal, sub- and super-diagonal matrices
for n_5 = 1 : n_layer
A=L(n_5).layer_npts;
if (n_5 == 1)
Main_diagonal1 = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2 = [ (L(n_5).H_centerI) , Main_diagonal1];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-1]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-1]); 
SubL1_diagonal2 = [ (L(n_5).H_minusI) , SubL1_diagonal1];
elseif (n_5 == n_layer)
Main_diagonal1 = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2 = [ (L(n_5).H_centerI) , Main_diagonal1];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-1]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-1]); 
SubL1_diagonal2 = [ (L(n_5).H_minusI) , SubL1_diagonal1];
elseif (L(n_5+1).layer_cntr == 1)    
Main_diagonal1  = repmat(L(n_5).H_centerJ,[1,A-2]);
Main_diagonal2  = [ (L(n_5).H_centerI) , Main_diagonal1 , (L(n_5+1).H_centerIL)];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-2]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1 , (L(n_5+1).H_minusI)' ];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-1]); 
SubL1_diagonal2 = [ (L(n_5).H_minusI) , SubL1_diagonal1];
elseif (L(n_5).layer_cntr == 1)
Main_diagonal1  = repmat(L(n_5).H_centerJ,[1,A-3]);
Main_diagonal2  = [(L(n_5).H_centerI),(L(n_5).H_centerIR), Main_diagonal1, L(n_5+1).H_centerIL ];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-2]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1 , (L(n_5+1).H_minusI)' ];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-2]);  
SubL1_diagonal2 = [ L(n_5).H_minusI , (L(n_5).H_plusI)' , SubL1_diagonal1];
elseif (L(n_5-1).layer_cntr == 1)
Main_diagonal1 = repmat(L(n_5).H_centerJ,[1,A-2]);
Main_diagonal2 = [(L(n_5).H_centerI) , L(n_5).H_centerIR , Main_diagonal1];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-1]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-2]); 
SubL1_diagonal2 = [ (L(n_5).H_minusI) , (L(n_5).H_plusI)' , SubL1_diagonal1]; 
else 
Main_diagonal1 = repmat(L(n_5).H_centerJ,[1,A-1]);
Main_diagonal2 = [ (L(n_5).H_centerI) , Main_diagonal1];
SupR1_diagonal1 = repmat(L(n_5).H_plusJ,[1,A-1]);
SupR1_diagonal2 = [ (L(n_5).H_plusI) , SupR1_diagonal1];
SubL1_diagonal1 = repmat(L(n_5).H_minusJ,[1,A-1]); 
SubL1_diagonal2 = [ (L(n_5).H_minusI) , SubL1_diagonal1];
end  
Main_diag{n_5}  = Main_diagonal2;
SupR1_diag{n_5} = SupR1_diagonal2;
SubL1_diag{n_5} = SubL1_diagonal2;   
end
Amd1  = cell2mat(Main_diag);
Amd2  = [Amd1  , infboundary_h1];
Asup1 = cell2mat(SupR1_diag);
Asup2 = [Asup1  , infboundary_h2];
Asup2(:,((8*tot_mesh_pts-7):8*tot_mesh_pts))=[];
Asub1 = cell2mat(SubL1_diag);
Asub1(:,9:16)=[];
Asub2 = [Asub1  , infboundary_h2];
Amd  = reshape(Amd2,[8,8,tot_mesh_pts+1]);
Asup = reshape(Asup2,[8,8,(tot_mesh_pts)]);
Asub = reshape(Asub2,[8,8,(tot_mesh_pts)]);
[p,q,n] = size(Amd); 
v = Amd(:);     % The main diagonal elements
if n>1          % then the sub and super diagonal blocks.
v=[v;Asub(:)];  % sub-diagonal
v=[v;Asup(:)];  % super-diagonal
end
% and now generate the index arrays. first the main diagonal
[ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-1);
rind = 1+ind1(:)+p*ind3(:);
cind = 1+ind2(:)+q*ind3(:);
% then the sub and super diagonal blocks;
% sub-diagonal
[ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-2);
rind = [rind;1+p+ind1(:)+p*ind3(:)];
cind = [cind;1+ind2(:)+q*ind3(:)];
% super-diagonal
rind = [rind;1+ind1(:)+p*ind3(:)];
cind = [cind;1+q+ind2(:)+q*ind3(:)];
% Finally build the final array all in one and call to sparse
Hmatrix = sparse(rind,cind,v,n*p,n*q);
setup_time=setup_time+toc; tic; % Counter timer for waitbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.5 Use the matlab function "eigs" to calculate eigenvalues and 
%%%% eigenfunctions for Conduction and Valence bands, and then order for
%%%% relavant subbands
%%%% Holes
[V_val_h,E_val_h] = eigs(Hmatrix ,n_Esub, Holes_Eval);
[Data(Cn_kpts).E_h,Data((Cn_kpts)).index_h] = sort((diag(E_val_h)),'descend');
V_val_h = V_val_h(:,Data((Cn_kpts)).index_h);
%%%% Electrons
[V_val_e,E_val_e] = eigs(Hmatrix ,n_Esub, Elect_Eval);
[Data(Cn_kpts).E_e,Data((Cn_kpts)).index_e] = sort((diag(E_val_e)),'ascend');
V_val_e = V_val_e(:,Data((Cn_kpts)).index_e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 2.5 Use the matlab function "eigs" to calculate eigenvalues and 
%%%% eigenfunctions for Conduction and Valence bands, and then order for
%%%% relavant subbands
for n_6=1:n_Esub
%%%% Seperate eigenfunctions into relevant states for electrons and holes
vectors_h      = reshape(V_val_h(:,n_6),8,tot_mesh_pts+1)';
V_ecu_h(:,n_6) = vectors_h(:,1);
V_ecd_h(:,n_6) = vectors_h(:,2);
V_hhu_h(:,n_6) = vectors_h(:,3);
V_hhd_h(:,n_6) = vectors_h(:,4);
V_lhu_h(:,n_6) = vectors_h(:,5);
V_lhd_h(:,n_6) = vectors_h(:,6);
V_sou_h(:,n_6) = vectors_h(:,7);
V_sod_h(:,n_6) = vectors_h(:,8);
vectors_e      = reshape(V_val_e(:,n_6),8,tot_mesh_pts+1)';
V_ecu_e(:,n_6) = vectors_e(:,1);
V_ecd_e(:,n_6) = vectors_e(:,2);
V_hhu_e(:,n_6) = vectors_e(:,3);
V_hhd_e(:,n_6) = vectors_e(:,4);
V_lhu_e(:,n_6) = vectors_e(:,5);
V_lhd_e(:,n_6) = vectors_e(:,6);
V_sou_e(:,n_6) = vectors_e(:,7);
V_sod_e(:,n_6) = vectors_e(:,8);
%%%% Calculate k=0 zone center eigenfunctions
if Cn_kpts == 1
V_ECu_h(:,n_6) = vectors_h(:,1).*conj(vectors_h(:,1));
V_ECd_h(:,n_6) = vectors_h(:,2).*conj(vectors_h(:,2));
V_HHu_h(:,n_6) = vectors_h(:,3).*conj(vectors_h(:,3));
V_HHd_h(:,n_6) = vectors_h(:,4).*conj(vectors_h(:,4));
V_LHu_h(:,n_6) = vectors_h(:,5).*conj(vectors_h(:,5));
V_LHd_h(:,n_6) = vectors_h(:,6).*conj(vectors_h(:,6));
V_SOu_h(:,n_6) = vectors_h(:,7).*conj(vectors_h(:,7));
V_SOd_h(:,n_6) = vectors_h(:,8).*conj(vectors_h(:,8));  
E_function_h(:,n_6) = V_ECu_h(:,n_6)+V_ECd_h(:,n_6)+V_HHu_h(:,n_6)+V_HHd_h(:,n_6)+V_LHu_h(:,n_6)+V_LHd_h(:,n_6)+V_SOu_h(:,n_6)+V_SOd_h(:,n_6);
V_ECu_e(:,n_6) = vectors_e(:,1).*conj(vectors_e(:,1));
V_ECd_e(:,n_6) = vectors_e(:,2).*conj(vectors_e(:,2));
V_HHu_e(:,n_6) = vectors_e(:,3).*conj(vectors_e(:,3));
V_HHd_e(:,n_6) = vectors_e(:,4).*conj(vectors_e(:,4));
V_LHu_e(:,n_6) = vectors_e(:,5).*conj(vectors_e(:,5));
V_LHd_e(:,n_6) = vectors_e(:,6).*conj(vectors_e(:,6));
V_SOu_e(:,n_6) = vectors_e(:,7).*conj(vectors_e(:,7));
V_SOd_e(:,n_6) = vectors_e(:,8).*conj(vectors_e(:,8));  
E_function_e(:,n_6) = V_ECu_e(:,n_6)+V_ECd_e(:,n_6)+V_HHu_e(:,n_6)+V_HHd_e(:,n_6)+V_LHu_e(:,n_6)+V_LHd_e(:,n_6)+V_SOu_e(:,n_6)+V_SOd_e(:,n_6);
end
end
Data(Cn_kpts).kx = kx;
Data(Cn_kpts).ky = ky;
eigs_time = eigs_time+toc;
waitbar((Cn_kpts/Tn_kpts),w);
end
end
end
close(w) 
disp('Total setup time'); disp(setup_time);
disp('Total eigs time'); disp(eigs_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 3.1 Check Dispersion relations
Cn_kpts = 0;
for y = 0:1:n_kpts-1
for x = 0:1:n_kpts-1 
if x >= y
Cn_kpts = Cn_kpts +1;  
%%%% Pull out CB and VB Dispersion energies in [01] and [11] directions
for w=1:n_Esub
if y == 0   
Eh_k_01(w,x+1)=Data(Cn_kpts).E_h(w);
Ee_k_01(w,x+1)=Data(Cn_kpts).E_e(w);
end
if x == y   
Eh_k_11(w,y+1)=Data(Cn_kpts).E_h(w);
Ee_k_11(w,y+1)=Data(Cn_kpts).E_e(w);
end
if x == 0 && y==0
Eh_subband(:,w) = ones(sum(layer_npts),1).*Eh_k_01(w,1)+E_function_h((1:sum(layer_npts)),w);
Ee_subband(:,w) = ones(sum(layer_npts),1).*Ee_k_01(w,1)+E_function_e((1:sum(layer_npts)),w);
end
end
end
end
end
k=0:(k_range/(n_kpts-1)):k_range;
kpar01= ((2*pi)/a0)*sqrt(k.*k);
kpar11= ((2*pi)/a0)*sqrt(k.*k + k.*k);
Z_profile = 1:sum(layer_npts);
E_electron  = [B_P(1).E_celectron B_P(2).E_celectron B_P(3).E_celectron B_P(4).E_celectron B_P(5).E_celectron ];
E_heavyhole = [B_P(1).E_heavyhole B_P(2).E_heavyhole B_P(3).E_heavyhole B_P(4).E_heavyhole B_P(5).E_heavyhole ];
E_lighthole = [B_P(1).E_lighthole B_P(2).E_lighthole B_P(3).E_lighthole B_P(4).E_lighthole B_P(5).E_lighthole ];
E_spinorbit = [B_P(1).E_spinorbit B_P(2).E_spinorbit B_P(3).E_spinorbit B_P(4).E_spinorbit B_P(5).E_spinorbit ];

E_electron_unp  = [B_P(1).E_celectron_unp B_P(2).E_celectron_unp B_P(3).E_celectron_unp B_P(4).E_celectron_unp B_P(5).E_celectron_unp ];
E_heavyhole_unp = [B_P(1).E_heavyhole_unp B_P(2).E_heavyhole_unp B_P(3).E_heavyhole_unp B_P(4).E_heavyhole_unp B_P(5).E_heavyhole_unp ];
E_lighthole_unp = [B_P(1).E_lighthole_unp B_P(2).E_lighthole_unp B_P(3).E_lighthole_unp B_P(4).E_lighthole_unp B_P(5).E_lighthole_unp ];
E_spinorbit_unp = [B_P(1).E_spinorbit_unp B_P(2).E_spinorbit_unp B_P(3).E_spinorbit_unp B_P(4).E_spinorbit_unp B_P(5).E_spinorbit_unp ];

Data_Strain_BO = [Z_profile' E_electron_unp' E_heavyhole_unp' E_lighthole_unp' E_spinorbit_unp' E_electron' E_heavyhole' E_lighthole' E_spinorbit'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 3.2 Plot (Figure 1) QW profile, the confining potential, 
%%%% subband levels, envelope functions and then (Figure 1) subband
%%%% dispersion in [01] and [11] directions
figure (1)
hold on
title('QW Subband energies and envelope functions')
grid on
plot(Z_profile, (E_electron' - Offset),'r', Z_profile, (E_heavyhole' - Offset),'r',...
     Z_profile, (E_lighthole' - Offset),'b', Z_profile, (E_spinorbit' - Offset),'g');
plot(Z_profile, (Eh_subband(:,:)'));
plot(Z_profile, (Ee_subband(:,:)'));
axis([150 350 -0.1 1.2]);
ylabel('Energy (eV)');
xlabel('z direction (nm)');  
hold off
figure (2)
hold on
title('Subband Dispersion')
grid on
plot(kpar01, (Eh_k_01)','r', -(kpar11), (Eh_k_11)','r');
plot(kpar01, (Ee_k_01)','r', -(kpar11), (Ee_k_11)','r');
axis([-0.06 0.06 -0.1 1.2]);
text(-0.06,-0.01,'[11]');
text(0.055,-0.01,'[01]');
ylabel('Energy (eV)');
xlabel('K|| (2 \pi /a_{0})');  
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stage 3.3 Define and Save data files
Data_QW_profile  = [Z_profile',(E_electron' - Offset ),(E_heavyhole' - Offset),(E_lighthole' - Offset),(E_spinorbit' - Offset)];
Data_QW_efunct   = [Z_profile',(Eh_subband(:,:)),(Ee_subband(:,:))]; 
Data_QW_dispersion=[kpar01',(Eh_k_01)',(Ee_k_01)',(kpar11)',(Eh_k_11)',(Ee_k_11)'];