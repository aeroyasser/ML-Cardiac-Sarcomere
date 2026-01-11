% -----------------------------------
% This code uses Rice model data to generate a shall/deep learning ML model
% using Artificial Neural Network
% ---------------------------------
% Yasser Aboelkassem, PhD
% University of Michigan - Flint
% September 24, 2025
% -------------------------------------------------------------

clc; clear all; close all;

%-----------------------------------
% Protocole # 1 (Dynamic Twitch Data)
%-----------------------------------
%-----------------------------------
% Step 1: Load Training Data 1 
%-----------------------------------
data_1=readmatrix('ML_Dynm_Rice_Protocol1_Fig5.csv');  % import data from your folder
x=data_1(1:5,:);  % input  Layer
y=data_1(6:8,:);  % output Layer
%------------------------
%----------------------------------
% Step 2: Normalize  Data if Needed
%----------------------------------

%-----------------------------------------
% Step 3: Artificial Neural Network Architecture
%-----------------------------------------
% You can optimize the number of hidden layers by calculating the rms as
% below
% let's start by a single hiddel layer

hiddenLayerSize=[10 10 5];         

net=fitnet(hiddenLayerSize);       % net object for ANN
net.divideParam.testRatio=70/100;  % use 70% on the input data as a training set OR split it ureself
net.divideParam.valRatio= 30/100;  % use 30% on the input data as a validation set
net.divideParam.testRatio= 0/100;  % use  0% on the input data as a test set
net.trainParam.epochs=1000;        % number of epoch
%-----------------------------------------
% Step 3: Train the ANN
%-----------------------------------------
% returns the trained network with correct weights and biases
% this can be used then to obtain the predicted output for a given new
% input set
[tarined_Net,tr]=train(net,x,y);  

%-----------------------------------------
% Step 4: Performance the ANN
%-----------------------------------------

yTrain     =tarined_Net(x(:,tr.trainInd));  % predicted output based on percentage trainned  data only
yTrainTrue =y(:,tr.trainInd);   % Actual output that we kept for validation


yValTrue=y(:,tr.valInd);   % Actual output that we kept for validation
yVal  =tarined_Net(x(:,tr.valInd));    % predicted output based on percentage validation data only

% Check root mean squared error
rmsTrain= sqrt(mean((yTrain-yTrainTrue).^2));  % calculate root mean squared error
rmsValid= sqrt(mean((yVal-yValTrue).^2));  % calculate root mean squared error
%-----------------------------
%-----------------------------------------
% Step 5: Test my ANN
%-----------------------------------------
%------------------------------------------------------------------------------
% Rice Validation
% Load Test#1 
%------------------------------------------------------------------------------
Data_Test1=readmatrix('ML_Dynm_Rice_Protocol1_Test1.csv');  % import data 
x_test1=Data_Test1(1:5,:);   % input Layer
y_test1=Data_Test1(6:8,:);   % output Layer
Pred1=tarined_Net(x_test1);  % predicted results (NB: this is all outputs states 
%------------------------------------------------------------------------

%----------------------------------------
% Rice: Validation
% Load Test#2  
%----------------------------------------
Data_Test2=readmatrix('ML_Dynm_Rice_Protocol1_Test2.csv');  % import data 
x_test2=Data_Test2(1:5,:);    % input  layer
y_test2=Data_Test2(6:8,:);    % output layer
Pred2=tarined_Net(x_test2);   % predicted results (NB: this is all outputs states 
