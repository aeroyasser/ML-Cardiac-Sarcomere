% Authors: Yasser Aboelkassem 
% Year  :  Summer 2025
%-----------------------------
% This function calculates the Calcium transient 
%--------------------------------------------------------------------------
% Input| [time in ms] 
%-------
% t:     [ms] simulation current time
%--------
% Output| [micro Mole] 
%---------
% Ca2 :   Calcium transient interpolation 
%--------------------------------------------------------------------------%
function Ca = Ca_Transient(t,Ca_amplitude)
tnew=rem(t,1000);
%-------------------------------------------           
% Formula  ( Rice et al., 2008 Biophys J.)
%-------------------------------------------           
t_start=0;         % [ms]
tau1=20;           % [ms];
tau2=110;          % [ms]
%Ca_amplitue=1.45;  % 1.45, 1.25, 1.15, 1.05, 0.95 0.85 [micro Mole]
Ca_diastolic=0.09; % [micro Mole]
A=(tau1/tau2);
Beta=(A)^(-1/(A-1))-(A)^(-1/(1-1/A));
B=(Ca_amplitude-Ca_diastolic)/Beta;

Ca=Ca_diastolic.*(tnew<=t_start)+(tnew>t_start).*(B.*(exp(-(tnew-t_start)/tau1)-exp(-(tnew-t_start)/tau2))+Ca_diastolic);
end


