%-------------------------------------------------------------------------------%
%                    University of Michigan                                   % 
%                   Engineering and Technology                                %                         
%          -----------------------------------------------                      %
%                     Aboelkassem Lab                                           %
%-------------------------------------------------------------------------------%
% This code is given to reproduce the work J. J. Rice et al 2008
% Reference: Biophyscial J. Vol 92, 2008
%--------------------------------------------------------------------------
%-------------                                                            % 
%  Author: Yasser Aboelkassem, PhD                                        %
%  January 2025                                                           %  
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                                                                         %
%                     (Dynamic Simulations)                               %
%                             TWITCHES                                   %
%                                                                         %
%           i.e., simulations at time-dependent Ca-concentaions.          %
%-------------------------------------------------------------------------%
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

%------------------
% Input parameters
%-----------------
Protocol=1;                   % 1 = Fig 5 , 2= Fig 7
Animal = 1;                   % ( = 1 ) for "rat" ,and  ( = 2)for "rabbit" 
Cell_type =2;                 % ( = 1 ) for "isolated cell", ( = 2) "trabeculae cell"
KSE = 50;                      %(norm Force/um) series elastic element
Temp = 22.5;                  % Temperature (C)
F_aftload_const =0.0;          % (0:0.1:0.8) will be used when simulating Isotonic


SLo=[1.8 1.9 2 2.1 2.2 2.3];                    % sarcomere length (um)
Ca_amplitude = [0.85 0.95 1.05 1.15 1.25 1.45]; % maximum Ca [uM]
Numb_Ca=length(Ca_amplitude);
Numb_SLo=length(SLo);

for i=5:5  %1:Numb_SLo

   for j=6:6%1:Numb_Ca % to simulate variouse levels of Ca

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
SL_int      = SLo(i); %  Sarcomere length                                   % Eq. (38)
intf_int    = 0;      %  Integrated force                                      % Eq. (39) 
X1_int=0;
% vector for all ICs
y0=[N_NoXB_int; P_NoXB_int; N_int; P_int; XBprer_int; XBpostr_int;  SL_int; xXBpostr_int; xXBprer_int; TRPNCaL_int; TRPNCaH_int; X1_int];

%------------------------------
% ODE Solver (Yasser's version)
%------------------------------
Tp=100000; % time in ms
%options = odeset('RelTol', 1e-012, 'AbsTol', [1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012 1e-012]);
%options = odeset('RelTol',1e-8);
[Time,Y]=ode15s(@(t,y) Dynamic_Rice_Yas_MF_Model(t,y, Animal, Cell_type, Protocol,Ca_amplitude(j),SLo(i),KSE,Temp,F_aftload_const),[0  Tp],y0);

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
%intf    = Y(:,12); %  Integrated force                                     % Eq. (39) 
X1=Y(:,12);
%----------------------------------------------------------------------------

[F_active,F_passive,F_total,SOVFThick]=Dynamic_Active_Passive_Force(Animal,Cell_type,SL,xXBpostr,XBpostr,xXBprer,XBprer);

%---------------
% Active Twitch
%----------------
figure(1)
plot(Time,F_active,'g')

    end % finish j
end % finish i

figure(2)
plot(Time,SL)

toc

   
