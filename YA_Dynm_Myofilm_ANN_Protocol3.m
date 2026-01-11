% -----------------------------------
% This code uses Rice model data to generate a shall/deep learning ML model
% using Artificial Neural Network
% Good source:  https://www.youtube.com/watch?v=xOzh6PMk21I
% ---------------------------------
% Yasser Aboelkassem, PhD
% University of Michigan - Flint
% September 24, 2025
% -------------------------------------------------------------
tic
clc; clear all; close all;

%-----------------------------------
% Protocole # 3 (Dynamic Twitch Data) (Fig6 A & B)
%-----------------------------------
%-----------------------------------
% Step 1: Load Training Data 1 
%-----------------------------------
data_1=readmatrix('ML_Dynm_Rice_Protocol3_Fig6.csv');  % import data [Time,Temp,Ca,SLo,KSE, Fact, Ftot SL]
x=data_1(1:5,:);  % input  [Time,Temp,Ca,SLo,KSE]
y=data_1(7:8,:);  % output [Fact, Ftot SL]

%------------------------
% Step 2: Visualize  Data
%-------------------------
% figure(1)
% plot(x(1,:),y(1,:),'g-');  % x1 = Time
% figure(1)
% plot(x(1,:),y(2,:),'g-');  % x1 = Time
%----------------------------------
% Step 3: Normalize  Data if Needed
%----------------------------------

%-----------------------------------------
% Step 3: Artificial Neural Network Architecture
%-----------------------------------------
% You can optimize the number of hidden layers by calculating the rms as
% below
% let's start by a single hiddel layer

hiddenLayerSize=[5];            % 10 neurons in  first HL & 5 in second HL
% hiddenLayerSize=[10 10];         % two hiddel layer with 10 neurons each 

net=fitnet(hiddenLayerSize);       % net object for ANN
net.divideParam.testRatio=70/100;  % use 70% on the input data as a training set OR split it ureself
net.divideParam.valRatio= 30/100;  % use 30% on the input data as a validation set
net.divideParam.testRatio= 0/100;  % use  0% on the input data as a test set
net.trainParam.epochs=1000;          % number of epoch
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
% Plot output from validation
%-----------------------------
% figure(2)
% plot(yValTrue,yVal,'o')
% figure(3)
% plot(rmsTrain,'*')


%-----------------------------------------
% Step 5: Test my ANN
%-----------------------------------------
%------------------------------------------------------------------------------
% Rice: Fig7AB Validation
% Load Test#1 Data @ Temp= 22.5, KSE=1 & vary Camax, SL=1.9 
%------------------------------------------------------------------------------
Data_Test1=readmatrix('ML_Dynm_Rice_Protocol3_Test1_Fig6AB.csv');  % import data [Time,Temp,Ca,SLo,KSE, Fact, Ftot SL]
x_test1=Data_Test1(1:5,:);   % input  [Time,Temp,Ca,SLo,KSE]
y_test1=Data_Test1(7:8,:);   % output [Fact, Ftot SL]
Pred1=tarined_Net(x_test1);  % predicted results (NB: this is all outputs states 

Time=x_test1(1,:);

for i=200:200:1200
figure(10)
plot(Time(i-200+1:1:i),y_test1(2,i-200+1:1:i),'-o'); % SL Shortening
hold on
plot(Time(i-200+1:1:i),Pred1(2,i-200+1:1:i),'g');

figure(11)
plot(Time(i-200+1:1:i),y_test1(1,i-200+1:1:i),'-o'); % Total Force
hold on
plot(Time(i-200+1:1:i),Pred1(1,i-200+1:1:i)-3*10^-3,'g');
end

%------------------------------------------------------------------------

toc


