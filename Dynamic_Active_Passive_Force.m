%--------------------------------------------------------------------------
%                  University of Michigan
%                        Aboelkassem Lab
%           -----------------------------
%  Copyright: Yasser Aboelkassem, January 2025
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This code is given to reproduce thw work J. J. Rice et al 2008
% Reference: Biophyscial J. Vol 92, 2008
%--------------------------------------------------------------------------
function [F_active,F_passive,F_total, SOVFThick]=Dynamic_Active_Passive_Force(Animal,Cell_type,SL,xXBpostr,XBpostr,xXBprer,XBprer)

%--------------
% Parameters
%--------------
%--------------------
% Sarcomere Geometry
%--------------------
len_thin  = 1.2;      % (um) thin filament length
len_thick = 1.65;     % (um) thick filament length
len_hbare = 0.1;      % (um) length of bare portion of thick filament
%-------------------------------
% Thin filament regulation and crossbridge cycling
%---------------------------------------------------

fapp    = (Animal==1)*500e-3+(Animal==2)*500e-3;      % (1/ms) XB on rate
gapp    = (Animal==1)*70e-3+(Animal==2)*70e-3;      % (1/ms) XB off rate
hf      = (Animal==1)*2000e-3+(Animal==2)*2000e-3;    % (1/ms) rate between pre-force and force states
hb      = (Animal==1)*400e-3+(Animal==2)*400e-3; % (1/ms) rate between pre-force and force states
gxb     = (Animal==1)*70e-3+(Animal==2)*70e-3;    %  (1/ms) ATP consuming transition rate : ( = 1) "rat", (=2) "rabbit" 


% fapp    = 500e-3;     % (1/ms) XB on rate
% gapp    = 70e-3;      % (1/ms) XB off rate
% hf      = 2000e-3;    % (1/ms) rate between pre-force and force states
% hb      = 400e-3;     % (1/ms) rate between pre-force and force states
% gxb     = 70e-3;      % (1/ms) ATP consuming transition rate
x_0     = 0.007;      % (um) strain induced by head rotation
kxb     = 120;        % (mN/mm^2) maximal force
%----------------------------------------------
%------------------------------------------
% Compute single overlap fractions "SOVF"
%------------------------------------------
sovr_ze   = min(len_thick/2,SL/2);               % z-line end                 % Eq.(42)
sovr_cle  = max(SL/2-(SL-len_thin),len_hbare/2); % centerline of end          % Eq.(43)
len_sovr  = sovr_ze-sovr_cle;                    % single overlap length      % Eq.(44)
SOVFThick = len_sovr*2/(len_thick-len_hbare);    % thick flmnt overlap frac   % Eq.(45)
%--------------------------------------------------------------------------
%  Active Force
%---------------
%SSXBprer = (hb*fapp+gxb*fapp)/(gxb*hf+fapp*hf+gxb*gapp+hb*fapp+hb*gapp+gxb*fapp);             % Eq.(33)
SSXBpostr = (fapp*hf)/(gxb*hf+fapp*hf+gxb*gapp+fapp*hb+gapp*hb+gxb*fapp);                       % Eq.(34)
% normalization for scaling active and passive force (maximal force)
F_active = SOVFThick.*(xXBpostr.*XBpostr+xXBprer.*XBprer)./(x_0*SSXBpostr);                                                          
%--------------------------------------------------------------------------
%---------------
% Passive Force
%---------------
%-------------------------------------
% Normalized active and passive force
%------------------------------------
SLrest  = 1.9;        %(um) rest SL length for 0 passive force
PCon_t  = 0.002;      % (norm Force) passive force due to titin
PExp_t  = 10;         %   these apply to trabeculae and single cells only
SL_c    = 2.25;       % (um) resting length for collagen
PCon_c  = 0.02;       % (norm Force) passive force due to collagen
PExp_c  = 70;         %  these apply to trabeculae and single cells only
F_passive_t = sign(SL-SLrest).*PCon_t.*(exp(PExp_t.*abs(SL-SLrest))-1);         % Titin force term       % Eq.(47)
F_passive_c = heav(SL-SL_c).*PCon_c.*(exp(PExp_c.*abs(SL-SL_c))-1);             % Collagen force term    % Eq.(48)
F_passive = (Cell_type==1)*F_passive_t + (Cell_type==2)*(F_passive_t+F_passive_c);                                % Total passive force    % Eq.(49)
%--------------------------------------------------------------------------
F_total=(F_active+F_passive);
return
function x = heav(val)
if val >= 0
    x = 1;
else
    x =0;
end
return

