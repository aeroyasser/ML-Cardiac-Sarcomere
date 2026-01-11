%-------------------------------------------------------------------------------%
%                 University of Michigan- Flint                                 % 
%                     Engineering & Technology                                  %                         
%          -----------------------------------------------                      %
%                     Aboelkassem's Lab                                         %
%-------------------------------------------------------------------------------%
% This code is given to reproduce the work J. J. Rice et al 2008
% Reference: Biophyscial J. Vol 92, 2008
%--------------------------------------------------------------------------
%-------------                                                            % 
%  Author: Yasser Aboelkassem, PhD                                        %
%  Septeber 2025                                                           %  
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                                                                         %
%                     (Dynamic Simulations)                               %
%                             TWITCHES                                   %
%                                                                         %
%           i.e., simulations at time-dependent Ca-concentaions.          %
%-------------------------------------------------------------------------%
%
% Protocol 1: to generate Fig5(A & B) from Rice Paper using ML
%
%-------------------------------------------------------
clc; clear all; close all;
tic

%--------------------------------------------------------------------------
%
%                          Parameters guide
%
%--------------------------------------------------------------------------
%-------------------------------------------------------------
% Animal Type   ( rat , rabbit )
%-------------------------------------------------------------
% (Animal = 1 ) for "rat" , ( = 2)for "rabbit"

%------------
% Cell Type
%------------
% (Cell_type =1) for "isolated cell", ( = 2) "trabeculae cell"

%--------------------------------------------------------
% Muscle Protocol (Isoscarcometric, Isometric, Isotonic)
%--------------------------------------------------------
% If (Protocol = 1) Isocarcometric (Fig.5A-B)
% If (Protocol = 2) Isometric with internal shortening (Fig.7A-B)
%--------------------------------------------------------------------------
%                    
%                        Study cases guide
%
%--------------------------------------------------------------------------
% Case 1:   Isosarcometric Twitches for different SLo =const 
%                     "Generate Fig-5A"
%            --------------------------         
%-----------
% (Fig. 5A)
%-----------
% 1- Protocol = 1
% 2- Set SLo =[1.8 1.9 2.0 2.1 2.2 2.3] [um]
% 3- Set Ca_amplitude=1.45 [uM]; 
%-----------
% (Fig. 5B)
%-----------
% 1- Protocol = 1
% 2- Set SLo =2.3 [um] 
% 3- Set Ca_amplitude=[0.85 0.95 1.05 1.15 1.25 1.45] [uM]
%--------------------------------------------------------------------------
%------------------------------------------------------------------------
% Use the following to generate data for the Machine Learning Algorithm
%------------------------------------------------------------------------
%------------------
% Output from ANN
%------------------
F_active_output=[];       % Prepare array to append F_active output vector
F_total_output=[];        % Prepare array to append F_total  output vector
SL_output=[];             % Prepare array to append SL input output vector
%--------------
% Input to ANN
%--------------
Time_input=[];                % Prepare array to append Time 
Temp_input=[];                % Prepare array to append Temperature input vector
Ca_input=[];                  % Prepare array to append Ca input vector
SLo_input=[];                 % Prepare array to append SL input vector
KSE_input=[];                   % Prepare array to append Stiffness input vector
%--------------------------------
T            = 22.5;%(20:1:30);      % 22.5 Temperature in (C)
Ca_Ampl      = 1.45;%[0.85 0.95 1.05 1.15 1.25 1.45];%(0.80:0.025:1.5); %[0.85 0.95 1.05 1.15 1.25 1.45]; % maximum Ca [uM]
SLo_values   = [1.8 1.9 2 2.1 2.2 2.3]; %(1.8:0.05:2.3); %2.2 % initial Sarcomere length (um)
K            = 1;%linspace(1,50,10);       % Stiffness series elastic (norm Force/um)      
F_AL         = 0.5;%(0:0.1:0.85);    % Force afterload be used when simulating Isotonic

%--------------------------
% Prepare some indices
%-------------------------
i_Temp_max =length(T);
i_K_max    =length(K);     
i_SLo_max  =length(SLo_values);
i_Ca_max   =length(Ca_Ampl); 

%------------------
% Input parameters
%-----------------
Protocol=1;                   % 1 = Fig 5 , 2= Fig 7, 3= Fig6
Animal = 1;                   % ( = 1 ) for "rat" ,and  ( = 2)for "rabbit" 
Cell_type =2;                 % ( = 1 ) for "isolated cell", ( = 2) "trabeculae cell"
F_aftload_const =0.5;          % (0:0.1:0.8) will be used when simulating Isotonic
%---------------------
% Start Loop Here
%---------------------
for M1 = 1:i_Temp_max        % loop over Temperature
    for M2=1:i_K_max          % loop over Stiffness [1-50]  
        for M3 = 1:i_SLo_max  % loop over SL 
        for M4=1:i_Ca_max      % to simulate variouse levels of Camax
                
                %-----------------
                % Initial Values
                %-----------------
                Temp = T(M1);                 % Temperature (C)
                KSE  = K(M2);                  %(norm Force/um) series elastic element
                SLo  = SLo_values(M3) ;        % initial sarcomere length (um)
                Ca_amplitude= Ca_Ampl(M4);    %  
                %-------------------
                
                
%-----------------------------
% Step1: Initial Concditions:
%-----------------------------
N_NoXB_int  = 0.99;   %  Nonpermissive state in the absence of nearby XBs      % Eq. (7)
P_NoXB_int  = 0.01;   %  Permissive state in the absence of nearby XBs         % Eq. (8)
N_int       = 0.97;   %  Nonpermissive state  i.e. prevents strong bound XBs   % Eq. (14)
P_int       = 0.01;   %  Permissive state  i.e., allows strong bound XBs       % Eq. (15)
XBprer_int  = 0.01;   %  XBs in   strong bound dtstr but no rotation yet       % Eq. (16)
XBpostr_int = 0.01;   %  XBs in strong bound state and allowed to rotate       % Eq. (17)
xXBpostr_int= 0.0;    %  averaged distortion of XBpostr                        % Eq. (30)             
xXBprer_int = 0.007;  %  averaged distortion of XBprer state                   % Eq. (29)
TRPNCaL_int = 0.0;    %  Low-affinity Ca-bound troponin (uM)                   % Eq. (2)
TRPNCaH_int = 0.0;    %  High-affinity Ca-bound troponin (uM)                  % Eq. (1)
SL_int      = SLo; %  Sarcomere length                                   % Eq. (38)
intf_int    = 0;      %  Integrated force                                      % Eq. (39) 
X1_int=0;
% vector for all ICs
y0=[N_NoXB_int; P_NoXB_int; N_int; P_int; XBprer_int; XBpostr_int;  SL_int; xXBpostr_int; xXBprer_int; TRPNCaL_int; TRPNCaH_int; X1_int];

%------------------------------
% ODE Solver (Yasser's version)
%------------------------------
Tp=30000; % Time: 30 cycles (each cycle = 1000 ms)
%options = odeset('RelTol', 1e-012, 'AbsTol', [1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012]);
%options = odeset('RelTol',1e-8);
[Time,Y]=ode15s(@(t,y) ML_Dynamic_Rice_Yas_MF_Model(t,y, Animal, Cell_type, Protocol,Ca_amplitude,SLo,KSE,Temp,F_aftload_const),[0  Tp],y0);

%----------------
% Data Output
%--------------
N_NoXB  = Y(:,1);   %  Nonpermissive state in the absence of nearby XBs      % Eq. (7)
P_NoXB  = Y(:,2);   %  Permissive state in the absence of nearby XBs         % Eq. (8)
N       = Y(:,3);   %  Nonpermissive state  i.e. prevents strong bound XBs   % Eq. (14)
P       = Y(:,4);   %  Permissive state  i.e., allows strong bound XBs       % Eq. (15)
XBprer  = Y(:,5);   %  XBs in   strong bound dtstr but no rotation yet       % Eq. (16)
XBpostr = Y(:,6);   %  XBs in strong bound state and allowed to rotate       % Eq. (17)
SL      = Y(:,7);   %  Sarcomere length                                      % Eq. (38)
xXBpostr= Y(:,8);   %  averaged distortion of XBpostr                        % Eq. (30)             
xXBprer = Y(:,9);   %  averaged distortion of XBprer state                   % Eq. (29)
TRPNCaL = Y(:,10);  %  Low-affinity Ca-bound troponin (uM)                   % Eq. (2)
TRPNCaH = Y(:,11);  %  High-affinity Ca-bound troponin (uM)                  % Eq. (1)
%intf    = Y(:,12); %  Integrated force                                      % Eq. (39) 
X1=Y(:,12);
%----------------------------------------------------------------------------

[F_active,F_passive,F_total,SOVFThick]=ML_Dynamic_Active_Passive_Force(Animal,Cell_type,SL,xXBpostr,XBpostr,xXBprer,XBprer);

%---------------------
% Active Twitch Force
%---------------------
N_intp=200;                        % Number of interpolated points
Time_intp=linspace(Tp-1000,Tp,N_intp); % Interpolation Time points
Factive_intp=interp1(Time,F_active,Time_intp); % Interpolation Avtive Force points

% figure(1)
% plot(Time,F_active,'g')
% hold on
% plot(Time_intp,Factive_intp,'r o')
%-------------------
% Total Twitch Force
%---------------------
Ftotal_intp=interp1(Time,F_total,Time_intp); % Interpolation Total Force points
% figure(2)
% plot(Time,F_total,'g')
% hold on
% plot(Time_intp,Ftotal_intp,'r o')
%---------------------
%-------------------
% Sarcommere Length
%---------------------
SL_intp=interp1(Time,SL,Time_intp); % Interpolation Total Force points
% figure(3)
% plot(Time,SL,'b')
% hold on
% plot(Time_intp,SL_intp,'r o')
%------------------
% Return Ca values
%------------------
T_Ca=Time_intp-(Tp-1000);
Ca_append = Ca_Transient(T_Ca,Ca_amplitude);
% figure(4)
% plot(T_Ca,Ca_append)
%-------------

%------------------------------------
% Prepare Data in Vectors for Export
%------------------------------------
Temp_append     = Temp.*ones(1,N_intp);
KSE_append      = KSE.*ones(1,N_intp);
SL_append       = SLo.*ones(1,N_intp);
%--------------------------
% Input as a single vector
%-------------------------- 
Time_input =[Time_input,(Time_intp-(Tp-1000))./1000];
Temp_input =[Temp_input, Temp_append];                  
KSE_input  =[KSE_input, KSE_append]; 
SLo_input  =[SLo_input, SL_append];
Ca_input  =[Ca_input, Ca_append];
%--------------------------
% Output as a single vector
%-------------------------- 
F_active_output  =[F_active_output, Factive_intp]; % Append F_active output vector
F_total_output   =[F_total_output, Ftotal_intp];  % Append F_total output vector
SL_output        =[SL_output, SL_intp];           % Append SL output vector

    end % finish M4
    end % finish M3
    end % finish M2
end % finish M1



% ---------------
% Force Active and Force Total Twitch and SL
% ----------------
figure(5)
plot(Time_input,F_active_output,'b')
figure(6)
plot(Time_input,F_total_output,'r')
figure(7)
plot(Time_input,SL_output,'r')
figure(8)
plot(Time_input,Ca_input,'r')
%-------------------

% %-------------------------
% % Export Training data Protocol # 1
% %-------------------------
% ML_Dynamic_Rice_Protocol1_Fig5=[Time_input; Temp_input; KSE_input; SLo_input; Ca_input; F_active_output; F_total_output; SL_output];
% csvwrite('ML_Dynm_Rice_Protocol1_Fig5.csv',ML_Dynamic_Rice_Protocol1_Fig5)
% %------------------------------

% %------------------------------------------------
% % Export Protocol1, Test#1 @ Temp= 22.5, KSE=1 & Camax = 1.45
% %------------------------------------------------
% ML_Dynamic_Rice_Protocol1_Test1_Fig5A=[Time_input; Temp_input; KSE_input; SLo_input; Ca_input; F_active_output; F_total_output; SL_output];
% csvwrite('ML_Dynm_Rice_Protocol1_Test1_Fig5A.csv',ML_Dynamic_Rice_Protocol1_Test1_Fig5A)
% %------------------------------

%------------------------------------------------
% Export Protocol1 Test#2 @ Temp= 22.5, KSE=1 & SL = 2.2 
%------------------------------------------------
ML_Dynamic_Rice_Protocol1_Test2_Fig5B=[Time_input; Temp_input; KSE_input; SLo_input; Ca_input; F_active_output; F_total_output; SL_output];
csvwrite('ML_Dynm_Rice_Protocol1_Test2_Fig5B.csv',ML_Dynamic_Rice_Protocol1_Test2_Fig5B)
%------------------------------
toc
   