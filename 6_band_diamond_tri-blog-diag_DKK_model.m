clc
clear all

clc
clear all
%%%%  Set up scaling constant
% Define fundamental constants in S.I units. All energies are in eV from
% CODATA 2006 at http://physics.nist.gov/cuu/Constants/index.html
hbar=1.054571628E-34;               %Planck's constant [Js]
m0=9.10938215E-31;                  %Free electron mass[Kg]
e=1.602176487E-19;                  %Electron charge [C]
L=1E-10;                            %Length Scale [Angstrom]
Const=(hbar*hbar)/(2*m0*e*L*L);     %Scaling constant

%%%% Input parameters
k_range=0.2;
n_kpts=300;
%%%% Define Material constants
%%%% Zone center energy seperations
Eg7p=-0.29;  
Eg8p=0;  
Eg6m=0.89; 
%%%% Material parameters              
a0  = 5.6579;
L1  =13.38;
L2  =4.24;
L2m =-4.24;
L3  =5.69;
M_e =(1/0.038);
M_so=(1/0.095);
Lb  =((1/2)*(L2+L3));
Mu  =((1/2)*(L3-L2));
Lbm  =((1/2)*(L2m+L3));
Mum  =((1/2)*(L3-L2m));
L2n =1;
L3n =8;
L2n =6.00;
L3n =8.05;
Lbn =((1/2)*(L2n+L3n));
Mun =((1/2)*(L3n-L2n));
E_p = 25.49;
P_0 =(1/Const)*sqrt((Const)*E_p);
%%%% Dresselhaus Parameters
L=31.8;
M=5.1;
N=32.1;
%%%% Direct Interactions
%ZD1_G6mG8p   =(1/Const)*sqrt((Const)*E_pg6mg8p)*1;

V_potS=diag([ Eg8p Eg8p Eg8p ]);   
V_potD=diag([ Eg8p Eg8p Eg8p Eg8p Eg8p-0.001  Eg8p-0.001 ]);
for n=1:4   
Cn_kpts=0;
for k_step = 0:1:n_kpts-1 
    
if n == 1
x = k_step;
y = 0;
z = 0;
elseif n == 2
x = k_step;
y = k_step;
z = k_step;
elseif n == 3
x = (k_range*(n_kpts-1));
y = k_step;
z = 0;
elseif n == 4
x = k_step;
y = k_step;
z = 0;
end    
    
Cn_kpts=Cn_kpts+1;
kx = x*(k_range/(n_kpts-1));
ky = y*(k_range/(n_kpts-1));
kz = z*(k_range/(n_kpts-1));
k   = (2*pi/a0)*sqrt(kx*kx + ky*ky + kz*kz);
k_p = (2*pi/a0)*(kx + 1i*ky);
k_m = (2*pi/a0)*(kx - 1i*ky);
k_x  = (2*pi/a0)*kx;
k_y  = (2*pi/a0)*ky;
k_z  = (2*pi/a0)*kz;
%%%% Single group Hamiltonian
HS    = [  (k*k) 0 0
            0 (k*k) 0
            0 0 (k*k)];
                 
  
H_DKK  = [(L*k_x*k_x + M*(k_y*k_y + k_z*k_z))  (N*k_x*k_y)                     (N*k_x*k_z)
         (N*k_x*k_y)                      (L*k_y*k_y + M*(k_x*k_x + k_z*k_z)) (N*k_y*k_z)                                   
         (N*k_x*k_z)                      (N*k_y*k_z)                     (L*k_z*k_z + M*(k_x*k_x + k_y*k_y))  ];
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Adapted double group
HD    = [  (k*k) 0 0 0 0 0
            0 (k*k) 0 0 0 0
            0 0 (k*k) 0 0 0
            0 0 0 (k*k) 0 0
            0 0 0 0 (k*k) 0 
            0 0 0 0 0 (k*k)];
        
H_DKK_double =  [(L*k_x^2)/6 + (L*k_y^2)/6 + (2*L*k_z^2)/3 + (5*M*k_x^2)/6 + (5*M*k_y^2)/6 + (M*k_z^2)/3,                                                   -(3^(1/2)*N*k_z*(k_x - k_y*1i))/3,  3^(1/2)*((L*k_x^2)/6 - (L*k_y^2)/6 - (M*k_x^2)/6 + (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                                                       0,                                                  -(2^(1/2)*N*k_z*(k_x + k_y*1i))/2,                                    -(2^(1/2)*(L - M)*(k_x^2 + k_y^2 - 2*k_z^2))/6
                                                        -(3^(1/2)*N*k_z*(k_x + k_y*1i))/3,                    (L*k_x^2)/2 + (L*k_y^2)/2 + (M*k_x^2)/2 + (M*k_y^2)/2 + M*k_z^2,                                                                                  0,       3^(1/2)*((L*k_x^2)/6 - (L*k_y^2)/6 - (M*k_x^2)/6 + (M*k_y^2)/6 + (N*k_x*k_y*1i)/3), 6^(1/2)*((L*k_x^2)/6 - (L*k_y^2)/6 - (M*k_x^2)/6 + (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                  -(6^(1/2)*N*k_z*(k_x + k_y*1i))/6
      -3^(1/2)*((L*k_y^2)/6 - (L*k_x^2)/6 + (M*k_x^2)/6 - (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                                                  0,                    (L*k_x^2)/2 + (L*k_y^2)/2 + (M*k_x^2)/2 + (M*k_y^2)/2 + M*k_z^2,                                                         (3^(1/2)*N*k_z*(k_x - k_y*1i))/3,                                                  -(6^(1/2)*N*k_z*(k_x - k_y*1i))/6, 6^(1/2)*((L*k_y^2)/6 - (L*k_x^2)/6 + (M*k_x^2)/6 - (M*k_y^2)/6 + (N*k_x*k_y*1i)/3)
                                                                                       0, -3^(1/2)*((L*k_y^2)/6 - (L*k_x^2)/6 + (M*k_x^2)/6 - (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                    (3^(1/2)*N*k_z*(k_x + k_y*1i))/3, (L*k_x^2)/6 + (L*k_y^2)/6 + (2*L*k_z^2)/3 + (5*M*k_x^2)/6 + (5*M*k_y^2)/6 + (M*k_z^2)/3,                                     (2^(1/2)*(L - M)*(k_x^2 + k_y^2 - 2*k_z^2))/6,                                                 (2^(1/2)*N*k_z*(k_y + k_x*1i)*1i)/2
                                                       (2^(1/2)*N*k_z*(k_y + k_x*1i)*1i)/2, -6^(1/2)*((L*k_y^2)/6 - (L*k_x^2)/6 + (M*k_x^2)/6 - (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                   -(6^(1/2)*N*k_z*(k_x + k_y*1i))/6,                                           (2^(1/2)*(L - M)*(k_x^2 + k_y^2 - 2*k_z^2))/6,                                             ((L + 2*M)*(k_x^2 + k_y^2 + k_z^2))/3,                                                                                 0
                                          -(2^(1/2)*(L - M)*(k_x^2 + k_y^2 - 2*k_z^2))/6,                                                   -(6^(1/2)*N*k_z*(k_x - k_y*1i))/6, -6^(1/2)*((L*k_x^2)/6 - (L*k_y^2)/6 - (M*k_x^2)/6 + (M*k_y^2)/6 + (N*k_x*k_y*1i)/3),                                                        -(2^(1/2)*N*k_z*(k_x + k_y*1i))/2,                                                                                 0,                                             ((L + 2*M)*(k_x^2 + k_y^2 + k_z^2))/3];
 
      
      
      
%  HS    = [ (HS_C)       (zeros(2,4))  (zeros(2,2))
%            (zeros(4,2)) (-HS_V)       (zeros(4,2))
%            (zeros(2,2)) (zeros(2,4))  (-HS_S)    ];
       
HSingle = Const.*-(HS+H_DKK) + V_potS;
[V_valS,E_valS] = eig(HSingle);
[E_valS,index_] = sort((diag(E_valS)),'descend');
V_valS          = V_valS(:,index_);
E_valsS(:,Cn_kpts,n)=E_valS;

HDouble = Const.*-(HD+H_DKK_double) + V_potD;
[V_valD,E_valD] = eig(HDouble);
[E_valD,index_] = sort((diag(E_valD)),'descend');
V_valD          = V_valD(:,index_);
E_valsD(:,Cn_kpts,n)=E_valD;
end   
Cn_kpts1=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0:(k_range/(n_kpts-1)):k_range;
kpar001_X= ((2*pi)/a0)*sqrt(k.*k);
kpar011_U= ((2*pi)/a0)*sqrt(k.*k+k_range*k_range);
kpar011_K= ((2*pi)/a0)*sqrt(k.*k+k.*k);
kpar111_L= ((2*pi)/a0)*sqrt(k.*k + k.*k + k.*k);
%%%% Plot figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (1)

min_E=-0.6;
max_E=0.2;
k_X=0.2;
k_L=-0.2;

hold on
%set(gcf,'Position',[200 400 800 400])
plot(kpar001_X, (E_valsS(:,:,1)),'k:','LineWidth',3);
plot( -kpar111_L, (E_valsS(:,:,2)),'k:','LineWidth',3);
 plot(kpar001_X, (E_valsD(:,:,1)),'k','LineWidth',3);
 plot( -kpar111_L, (E_valsD(:,:,2)),'k','LineWidth',3);
 axis([k_L k_X min_E max_E]);
 %box(axes1,'on');
 xlabel('k a_0 / 2\pi');    %  label the x-axis
 ylabel('Energy (eV)');  %  label the y-axis
plot(zeros(1,10),linspace(min_E,max_E,10),'k')

% Create textarrow
annotation(figure(1),'textarrow',[0.207808564231738 0.163727959697733],...
    [0.157824933687003 0.156498673740053],'TextEdgeColor','none',...
    'TextLineWidth',1,...
    'FontName','Times New Roman',...
    'String',{'L[111]'},...
    'LineWidth',1);


% Create textarrow
annotation(figure(1),'textarrow',[0.82367758186398 0.874055415617128],...
    [0.155126482213439 0.153846153846154],'TextEdgeColor','none',...
    'TextLineWidth',1,...
    'FontName','Times New Roman',...
    'String',{'X[100]'},...
    'LineWidth',1);

% Create textbox
annotation(figure(1),'textbox',...
    [0.526188916876574 0.139257294429708 0.0279672544080605 0.0331564986737398],...
    'String',{'\Gamma'},...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[0.972549019607843 0.972549019607843 0.972549019607843]);
%CREATEAXES(PARENT1,X1,YMATRIX1,X2,YMATRIX2,X3,YMATRIX3)
%  PARENT1:  axes parent
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data
%  X3:  vector of x data
%  YMATRIX3:  matrix of y data

% %  Auto-generated by MATLAB on 25-May-2011 13:24:59
% 
% % Create axes
% axes1 = axes('Parent',Parent1);
% % Uncomment the following line to preserve the X-limits of the axes
% % xlim(axes1,[-0.2 0.2]);
% % Uncomment the following line to preserve the Y-limits of the axes
% % ylim(axes1,[-0.4 1.2]);
% box(axes1,'on');
% hold(axes1,'all');
% 
% % Create multiple lines using matrix input to plot
% plot1 = plot(X1,YMatrix1,'Color',[0 0 0],'Parent',axes1);
% set(plot1(1),'LineWidth',2);
% set(plot1(2),'LineWidth',2);
% set(plot1(3),'LineWidth',2);
% set(plot1(4),'LineWidth',2);
% set(plot1(5),'LineWidth',3);
% set(plot1(6),'LineWidth',3);
% set(plot1(7),'LineWidth',3);
% set(plot1(8),'LineWidth',3);
% 
% % Create multiple lines using matrix input to plot
% plot2 = plot(X2,YMatrix2,'Color',[0 0 0],'Parent',axes1);
% set(plot2(1),'LineWidth',2);
% set(plot2(2),'LineWidth',2);
% set(plot2(3),'LineWidth',2);
% set(plot2(4),'LineWidth',2);
% set(plot2(5),'LineWidth',3);
% set(plot2(6),'LineWidth',3);
% set(plot2(7),'LineWidth',3);
% set(plot2(8),'LineWidth',3);
% 
% % Create xlabel
% xlabel('k a_0 / 2\pi');
% 
% % Create ylabel
% ylabel('Energy (eV)');
% 
% % Create multiple lines using matrix input to plot
% plot(X3,YMatrix3,'Color',[0 0 0]);

%  axis([k_L k_X min_E max_E]);
%  %box(axes1,'on');
%  xlabel('k a_0 / 2\pi');    %  label the x-axis
%  ylabel('Energy (eV)');  %  label the y-axis
% plot(zeros(1,10),linspace(min_E,max_E,10),'k')


% for n=1:1:100
%     x=(n-1)/100;
% L1  =13.38 - 9.16*x;
% %L2  =-(4.24 - 3.86*x);
% L2  =-(4.24 - 4.62*x);
% m_hh = 1/(L1 + 2*L2);
% m_lh = 1/(L1 - 2*L2);
% R(n)=(m_hh/m_lh);
% Comp(n)=x;
% end
% figure(3)
% hold on
% plot(Comp,R)

% figure(2);
% clf
% set(gcf,'Position',[237 200 700 420])
% ylabel('Energy [eV]')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k111_L=fliplr(kpar111_L);
% k001_X=kpar001_X;
% k011_U=kpar011_U;
% k011_K=fliplr(kpar011_K);
% kdata=[k111_L k001_X k011_U k011_K];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %this mumbo-jumbo tranformation is necessary for plotting
% kplot=zeros(1,(4*n_kpts));
% kdiff=abs(diff(kdata));
% for i=2:4*n_kpts,
%   kplot(i)=sum(kdiff(1:i-1));
% end
% 
% EvalS_L_Gamma = (fliplr(E_valsS(:,:,2)));
% EvalS_Gamma_X = (E_valsS(:,:,1));
% EvalS_X_K     = (E_valsS(:,:,3));
% EvalS_K_Gamma = (fliplr(E_valsS(:,:,4)));
% 
% % EvalD_L_Gamma = (fliplr(E_valsD(:,:,2)));
% % EvalD_Gamma_X = (E_valsD(:,:,1));
% % EvalD_X_K     = (E_valsD(:,:,3));
% % EvalD_K_Gamma = (fliplr(E_valsD(:,:,4)));
% 
% Eval_Single = [EvalS_L_Gamma EvalS_Gamma_X EvalS_X_K EvalS_K_Gamma];
% %Eval_Double = [EvalD_L_Gamma EvalD_Gamma_X EvalD_X_K EvalD_K_Gamma];
% 
% maxE=max(max(Eval_Single));
% minE=min(min(Eval_Single));
%     
% hold on
% h=plot(kplot,Eval_Single','b'); 
% set(h,'MarkerSize',4)
% set(h,'MarkerFaceColor','blue')
% % h=plot(kplot,Eval_Double','r'); 
% % set(h,'MarkerSize',4)
% % set(h,'MarkerFaceColor','red')
% 
% 
% %draw vertical lines at major breakpoints
% 
% plot(kplot(n_kpts)*ones(1,10),linspace(minE,maxE,10),':')
% plot(kplot(2*n_kpts)*ones(1,10),linspace(minE,maxE,10),':')
% plot(kplot(3*n_kpts)*ones(1,10),linspace(minE,maxE,10),':')
% plot(kplot(4*n_kpts)*ones(1,10),linspace(minE,maxE,10),':')
% 
% [ax,h1,h2]=plotyy(kplot,zeros(4*n_kpts),kplot,zeros(4*n_kpts));
% hold off
% 
% axes(ax(1))
% axis([0 max(kplot) -30 10])
% axes(ax(2))
% axis([0 max(kplot) minE maxE])
% klabels=[k111_L k001_X k011_K];
% set(ax(1),'XTick',kplot)
% 
% set(ax(1),'XTickLabel',klabels)
% set(ax(1),'YTick',minE:10:maxE)
% set(ax(1),'XColor',[0 0 0])
% set(ax(1),'Box','off')
% set(ax(1),'YColor',[0 0 0])
% 
% set(ax(2),'XTick',kplot)
% set(ax(2),'XColor',[0 0 0])
% set(ax(2),'XAxisLocation','top')
% set(ax(2),'XTickLabel',[])
% set(ax(2),'YTick',0:10:maxe)
% set(ax(2),'YColor',[0 0 0])
% ylabel('energy [meV]')
% 
% write the high-symmetry point labels
% h=text(-.03,-5,'\Gamma');
% set(h,'FontSize',14)
% h=text(.97,-5,'X');
% set(h,'FontSize',14)
% h=text(1.33,-5,'K');
% set(h,'FontSize',14)
% h=text(2.38,-5,'\Gamma');
% set(h,'FontSize',14)
% 
% text(.42,-5,'[100]')
% text(1.78,-5,'[110]')
% text(2.8,-5,'[111]')


% 
% %%%%  Set up scaling constant
% % Define fundamental constants in S.I units. All energies are in eV from
% % CODATA 2006 at http://physics.nist.gov/cuu/Constants/index.html
% hbar=1.054571628E-34;               %Planck's constant [Js]
% m0=9.10938215E-31;                  %Free electron mass[Kg]
% e=1.602176487E-19;                  %Electron charge [C]
% L=1E-10;                            %Length Scale [Angstrom]
% Const=(hbar*hbar)/(2*m0*e*L*L);     %Scaling constant
% a0  = 5.6579;
% L1  =13.38;
% L2  =-4.24;
% L3  =5.69;
% Lb  =((1/2)*(L2+L3));
% Mu  =((1/2)*(L3-L2));
% k_range=0.1;
% n_kpts=300;
% % for z = 0:1:n_kpts-1
% c_kpts = 0;
% for y = 0:1:n_kpts-1
% for x = 0:1:n_kpts-1     
% if (x>=y)
% z=0;
% c_kpts=c_kpts+1;
% k_x = x*(k_range/(n_kpts-1));
% k_y = y*(k_range/(n_kpts-1));
% k_z = z*(k_range/(n_kpts-1));
% 
% k   = (2*pi/a0)*sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
% k_p = (2*pi/a0)*(k_x + 1i*k_y);
% k_m = (2*pi/a0)*(k_x - 1i*k_y);
% kx  = (2*pi/a0)*k_x;
% ky  = (2*pi/a0)*k_y;
% kz  = (2*pi/a0)*k_z;
% %V_potS=diag([0.1 -0.1 -0.1 0.1]); 
% V_potS=diag([-0.1 0 0 -0.1]);
% 
% HS_V  = [(k*k*(L1+L2)-kz*kz*(3*L2))            (-2*sqrt(3)*k_m*kz*L3)                (sqrt(3)*(-Lb*k_m*k_m + Mu*k_p*k_p))   (0)
%          (-2*sqrt(3)*k_p*kz*L3)                (k*k*(L1-L2)+kz*kz*(3*L2))            (0)                                    (-sqrt(3)*(-Lb*k_m*k_m + Mu*k_p*k_p))
%          (sqrt(3)*(-Lb*k_p*k_p + Mu*k_m*k_m))  (0)                                   (k*k*(L1-L2)+kz*kz*(3*L2))             (-2*sqrt(3)*k_m*kz*L3)
%          (0)                                   (-sqrt(3)*(-Lb*k_p*k_p + Mu*k_m*k_m)) (-2*sqrt(3)*k_p*kz*L3)                 (k*k*(L1+L2)-kz*kz*(3*L2)) ];
% 
% HS    = (-HS_V);  
%        
% HSingle = Const.*(HS) + V_potS;
% E_val = eig(HSingle);
% 
% HH(c_kpts) = E_val(4);
% LH(c_kpts) = E_val(1);
%  
% end
% end
% end
% 
% 
% Cn_kpts =0;
% for y = 0:1:n_kpts-1
% for x = 0:1:n_kpts-1 
%     if x >= y
%     Cn_kpts = Cn_kpts +1;
% 
%  %%% Calculate Joint Density of states for TE and TM polarisation      
% 
%         JDOS_TE((y+1),(x+1)) = HH(Cn_kpts);
%         JDOS_TM((y+1),(x+1)) = LH(Cn_kpts);
%        
%     end
% end
% end    
% 
%   JDOS_TE1 = ((JDOS_TE(:,:))+triu((JDOS_TE(:,:)),1)'); 
%   JDOS_TM1 = ((JDOS_TM(:,:))+triu((JDOS_TM(:,:)),1)');
% 
%   JDOS_TE2 = fliplr(JDOS_TE1);
%   JDOS_TM2 = fliplr(JDOS_TM1);
% 
%   JDOS_TE3 = horzcat(JDOS_TE2(:,1:(size(JDOS_TE2,2)-1)),JDOS_TE1);
%   JDOS_TM3 = horzcat(JDOS_TM2(:,1:(size(JDOS_TM2,2)-1)),JDOS_TM1);
%  
%   JDOS_TE4 = flipud(JDOS_TE3); 
%   JDOS_TM4 = flipud(JDOS_TM3);
% 
%  HH_Econtour(:,:) = vertcat(JDOS_TE4(1:(size(JDOS_TE4 ,1)-1),:),JDOS_TE3); 
%  LH_Econtour(:,:) = vertcat(JDOS_TM4(1:(size(JDOS_TM4 ,1)-1),:),JDOS_TM3);
%  
% figure(4)
% contour(HH_Econtour(1:599,1:599),'DisplayName','HH_Econtour(1:599,1:599)');figure(gcf)
% 
% figure(5)
% contour(LH_Econtour(1:599,1:599),'DisplayName','LH_Econtour(1:599,1:599)');figure(gcf)
%  
% figure(6)
% hold
% plot(HH(1:300),'r')
% plot(LH(1:300),'b')
 
 
 