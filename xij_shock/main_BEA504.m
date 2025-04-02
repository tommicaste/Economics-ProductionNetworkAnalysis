%% Compute result in "Sectoral Effect of Social Distancing" 
%(Table 2 and Figure 2)

%% Copyright: Jean-Noel Barrot, Basile Grassi and Julien Sauvagnat
%%AEA Paper and Procedings May 2021
%%Solving the model in BGS (AEA P&P) and performing exercice for various labor supply shocks
%%This version: the AEA P&P
%%Corresponding author: Basile Grassi

%% House keeping

clear all;
close all;


%% Load the library
addpath('lib')


%% Choose if you want the figure to pop-up
%if you want the figure to pop-up
set(0,'DefaultFigureVisible','on')
%if you don't want figure to pop-up
%set(0,'DefaultFigureVisible','off')


%% Calibrate

calibrate_BEA504_hat_fun


%% Exercice with Acatual Shocks

%%Choose the excel file name with the shocks
    excelfileshocks='Data_US_BEA405_DEC2020';


%%Baseline Calibration
    sigma=0.9;
    theta=0.5;
    epsilon=0.001;
    disp(table(sigma, theta,  epsilon));

    %computing the effect of shocks
    exercice_hat_fun_nicefigure(theta,sigma,epsilon,0.95,'US_BEA504_hat',excelfileshocks,0)

%%Cobb-Douglas
    sigma=1;
    theta=1;
    epsilon=1;
    disp(table(sigma, theta,  epsilon ));

    %computing the effect of shocks
    exercice_hat_fun_nicefigure(theta,sigma,epsilon,0.95,'US_BEA504_hat',excelfileshocks,0)
