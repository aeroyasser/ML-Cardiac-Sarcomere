% Implements myofilament model by Rice et al 2008

function dydt = Dynamic_Rice_Yas_MF_Model(t,y,Animal, Cell_type, Protocol,Ca_amplitude,SLo,KSE,Temp,F_aftload_const)

%-------------------
% Cal Ca2+
%-----------------
Ca = Ca_Transient(t,Ca_amplitude);

%-------------------
% Model Parameters  %
%-------------------
% Sarcomere Geometry
SLmax   = 2.4;        % (um) maximum sarcomere length
SLmin   = 1.4;        % (um) minimum sarcomere length
len_thin  = 1.2;      % (um) thin filament length
len_thick = 1.65;     % (um) thick filament length
len_hbare = 0.1;      % (um) length of bare portion of thick filament

% Temperature Dependence
Qkon = 1.5;
Qkoff = 1.3;
Qkn_p = 1.6;
Qkp_n = 1.6;
Qfapp = 6.25;
Qgapp = 2.5;
Qhf = 6.25;
Qhb = 6.25;
Qgxb = 6.25;

% Ca binding to troponin
kon     = 50e-3;      % (1/[ms uM])
koffL   = 250e-3;     % (1/ms)
koffH   = 25e-3;      % (1/ms)
perm50  = 0.5;        % perm variable that controls n to p transition
nperm   = 15;         %   in Hill-like fashion
kn_p    = (Animal==1)*500e-3+(Animal==2)*500e-3;     % (1/ms)
kp_n    = (Animal==1)*50e-3+(Animal==2)*50e-3;      % (1/ms)
% Thin filament regulation and crossbridge cycling
fapp    = (Animal==1)*500e-3+(Animal==2)*500e-3;      % (1/ms) XB on rate
gapp    = (Animal==1)*70e-3+(Animal==2)*70e-3;      % (1/ms) XB off rate
gslmod  = 6;          % controls SL effect on gapp
hf      = (Animal==1)*2000e-3+(Animal==2)*2000e-3;   % (1/ms) rate between pre-force and force states
hfmdc   = 5;          % 
hb      = (Animal==1)*400e-3+(Animal==2)*400e-3;     % (1/ms) rate between pre-force and force states
hbmdc   = 0;          % 
gxb     = (Animal==1)*70e-3+(Animal==2)*70e-3;    %  (1/ms) ATP consuming transition rate : ( = 1) "rat", (=2) "rabbit",  
sigmap  = 8;          % distortion dependence of STP using transition gxb
sigman  = 1;          % 
xbmodsp = (Animal==1)*1+(Animal==2)*4/3;    % mouse specific modification for XB cycling rates: ( = 1) "rat", (=2) "rabbit"
koffmod=1;
% Mean strain of strongly-bound states
x_0     = 0.007;      % (um) strain induced by head rotation
xPsi    = 2;          % scaling factor balancing SL motion and XB cycling

% Normalized active and passive force
SLrest  = 1.9 ;       % (um) rest SL length for 0 passive force
PCon_t  = 0.002;      % (norm Force) passive force due to titin
PExp_t  = 10;         %   these apply to trabeculae and single cells only
SL_c    = 2.25;       % (um) resting length for collagen
PCon_c  = 0.02;       % (norm Force) passive force due to collagen
PExp_c  = 70;         %   these apply to trabeculae and single cells only

% Calculation of complete muscle response
massf   = 0.00005e6;  % ([norm Force ms^2]/um) muscle mass
visc    = 0.003e3;    % ([norm Force ms]/um) muscle viscosity
kxb     = 120;        % (mN/mm^2) maximal force
Trop_conc = 70;       % (uM) troponin concentration

%---------------------
%     Variables      
%---------------------

% State Variables
N_NoXB  = y(1);   %  Nonpermissive state in the absence of nearby XBs    (Yasser) Eq. (7)
P_NoXB  = y(2);   %  Permissive state in the absence of nearby XBs       (Yasser) Eq. (8)
N       = y(3);   %  Nonpermissive state  i.e. prevents strong bound XBs (Yasser) Eq. (14)
P       = y(4);   %  Permissive state  i.e., allows strong bound XBs     (Yasser) Eq. (15)
XBprer  = y(5);   %  XBs in   strong bound dtstr but no rotation yet     (Yasser) Eq. (16)
XBpostr = y(6);   %  XBs in strong bound state and allowed to rotate     (Yasser) Eq. (17)
%P = 1-N-XBprer-XBpostr;       %  instead of using P as a state: Since, dp/dt = -(dN/dt + dXBprer/dt +dXBpostr/dt) 

SL      = y(7);   % Sarcomere length                           (Yasser) Eq. (38)
xXBpostr= y(8);   % averaged distortion of XBpostr             (Yasser) Eq. (30)             
xXBprer = y(9);   % averaged distortion of XBprer state        (Yasser) Eq. (29)
TRPNCaL = y(10);  % Low-affinity Ca-bound troponin (uM)        (Yasser) Eq. (2)
TRPNCaH = y(11);  % High-affinity Ca-bound troponin (uM)       (Yasser) Eq. (1)
X1      = y(12);  % integrated force
%--------------------------------------------------------------------------------------------------------
%---------------------
%     Equations      %
%---------------------

% Compute single overlap fractions
sovr_ze   = min(len_thick/2,SL/2);            % z-line end
sovr_cle=max(SL/2-(SL-len_thin),len_hbare/2); % centerline of end
len_sovr  = sovr_ze-sovr_cle;                 % single overlap length
SOVFThick = len_sovr*2/(len_thick-len_hbare); % thick filament overlap frac
SOVFThin  = len_sovr/len_thin;                % thin filament overlap frac

% Compute combined Ca binding to high- (w/XB) and low- (no XB) sites
Tropreg = (1-SOVFThin)*TRPNCaL + SOVFThin*TRPNCaH;
permtot = sqrt(1/(1+(perm50/Tropreg)^nperm));
inprmt = min(1/permtot, 100);

% Adjustments for Ca activation, temperature, SL, stress and strain
konT    = kon*Qkon^((Temp-37)/10);
koffLT  = koffL*Qkoff^((Temp-37)/10)*koffmod;
koffHT  = koffH*Qkoff^((Temp-37)/10)*koffmod;
kn_pT   = kn_p*permtot*Qkn_p^((Temp-37)/10);
kp_nT   = kp_n*inprmt*Qkp_n^((Temp-37)/10);
fappT   = fapp*xbmodsp*Qfapp^((Temp-37)/10);
gapslmd = 1 + (1-SOVFThick)*gslmod;
gappT   = gapp*gapslmd*xbmodsp*Qgapp^((Temp-37)/10);
hfmd    = exp(-sign(xXBprer)*hfmdc*((xXBprer/x_0)^2));
hbmd    = exp(sign((xXBpostr-x_0))*hbmdc*(((xXBpostr-x_0)/x_0)^2));
hfT     = hf*hfmd*xbmodsp*Qhf^((Temp-37)/10);
hbT     = hb*hbmd*xbmodsp*Qhb^((Temp-37)/10);
gxbmd   = heav(x_0-xXBpostr)*exp(sigmap*((x_0-xXBpostr)/x_0)^2)+...
  (1-heav(x_0-xXBpostr))*exp(sigman*(((xXBpostr-x_0)/x_0)^2));
gxbT    = gxb*gxbmd*xbmodsp*Qgxb^((Temp-37)/10);

% Regulation and corssbridge cycling state derivatives
dTRPNCaL  = konT*Ca*(1-TRPNCaL) - koffLT*TRPNCaL;
dTRPNCaH  = konT*Ca*(1-TRPNCaH) - koffHT*TRPNCaH;
dN_NoXB   = -kn_pT*N_NoXB + kp_nT*P_NoXB;
dP_NoXB   = -kp_nT*P_NoXB + kn_pT*N_NoXB;
dN        = -kn_pT*N + kp_nT*P;
% dP      = -kp_nT*P + kn_pT*N - fappT*P + gappT*XBprer + gxbT*XBpostr;
dXBprer   = fappT*P - gappT*XBprer - hfT*XBprer + hbT*XBpostr;
dXBpostr  = hfT*XBprer - hbT*XBpostr - gxbT*XBpostr;
dP        = -(dN+dXBprer+dXBpostr);

%--------------
% Active Force
%-------------
% steady-state fractions in XBprer and XBpostr using King-Altman rule
SSXBprer = (hb*fapp+gxb*fapp)/...
  (gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
SSXBpostr = fapp*hf/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);
% normalization for scaling active and passive force (maximal force)
Fnordv = kxb*x_0*SSXBpostr;
Non_norm_force = kxb*SOVFThick*(xXBpostr*XBpostr+xXBprer*XBprer);
F_active = Non_norm_force/Fnordv;

%---------------
% Passive Force
%---------------
ppforce_t = sign(SL-SLrest)*PCon_t*(exp(PExp_t*abs(SL-SLrest))-1); % titin
ppforce_c = heav(SL-SL_c)*PCon_c*(exp(PExp_c*abs(SL-SL_c))-1);     % collagen
F_passive = (Cell_type==1)*ppforce_t + (Cell_type==2)*(ppforce_t + ppforce_c);

%---------------
% Preload Force
%---------------
F_preload = sign(SLo-SLrest)*PCon_t*(exp(PExp_t*abs(SLo-SLrest))-1);

%-----------------
% Afterload Force
%-----------------
if     (Protocol==2)   % Isometric + elastic series (i.e., allow for shortening) (Fig3B)
       F_afterload = KSE*(SLo-SL);      
elseif (Protocol==4)   % Isotonic (Isosarcometric + fast rlease under const afterload force (Fig4A)
       F_afterload =F_aftload_const;
else
       F_afterload =0; %(Fig3C)    
end

%-------------
% Total Force
%-------------
F_total = (F_preload+F_afterload)-(F_passive+F_active);     % total force

% figure(11)
% plot(t,F_preload,'*r',t,F_afterload,'*b',t,F_passive,'*g',t,F_active,'*k')
% hold on
% %-------------
% figure(12)
% plot(t,F_total,'o')
% hold on


%-------------------------
%  Change in SL Equations
%  Yasser-SL 2-equations
%-------------------------

if (Protocol==1) || ((Protocol==4) && (t<=650)) % Isocarcometric OR Isotonic (constatnt length followed by fast release) Note: t is in [ms]
   dSL = 0;  % no change is sarcomere length
   dX1 = 0;
   
else
   dSL = X1; % isotonic case
   dX1 = ((F_total+(-X1)*visc)/massf)*...     % change in SL
      heav(SL-SLmin)*heav(SLmax-SL);
end
  
%--------------------------------
  
  
% Mean strain of strongly-bound states due to SL motion and XB cycling
dutyprer  = (hbT*fappT+gxbT*fappT)/...    % duty fractions using the
  (fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
dutypostr = fappT*hfT/...                 % King-Alman Rule    
  (fappT*hfT+gxbT*hfT+gxbT*gappT+hbT*fappT+hbT*gappT+gxbT*fappT);
dxXBprer = dSL/2+xPsi/dutyprer*(-xXBprer*fappT+(xXBpostr-x_0-xXBprer)*hbT);
dxXBpostr = dSL/2+xPsi/dutypostr*(x_0+xXBprer-xXBpostr)*hfT;

%-----------------------------------------------------
% Update state derivatives
dydt = zeros(size(y));
dydt(1)   = dN_NoXB;
dydt(2)   = dP_NoXB;
dydt(3)   = dN;
dydt(4)   = dP;
dydt(5)   = dXBprer;
dydt(6)   = dXBpostr;
dydt(7)   = dSL;
dydt(8)   = dxXBpostr;
dydt(9)   = dxXBprer;
dydt(10)  = dTRPNCaL;
dydt(11)  = dTRPNCaH;
dydt(12)  = dX1;
return
function x = heav(val)
if val >= 0
    x = 1;
else
    x =0;
end
return
