%%%%%%%%%%%%%% READ ME FILE 
%%%%%%%%%%%%%% MATLAB CODE OF "SECTORAL EFFECT OF SOCIAL DISTANCING" (Barrot, Grassi and, Sauvagnat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Folder structure

%lib contains some matlab scripts
%data_deep contains the raw data used by the matlab script
%fig contains the figures exported by the matlab script
%export contains the excel file produced by the matlab script
%data_matlab contains (mostly temporary) data produced by the matlab script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Raw data desciption

These files are in the data_deep folder. 

%data_deep/IOUse_Before_Redefinitions_PRO_DET.xlsx
	This files comes from the BEA and gives the input-output data.

%data_deep/Data_US_BEA405_DEC2020_laborsupply.xlsx
	This files contains the labor_supply shocks. 

	column A: sector code list
	column B: share_closed correspond to close_i in our paper
	column C: share_kid_constrained correspond to kids_i in our paper is the labor shock associated to our scenario (1)
	column D: share_nocompute correspond to wfh_i in our paper (or 1-wfh_i)
	column E: share_closednocomputerkid is the labor shock associated to our baseline scenario (3)
	column F: share_closednocomputer is the labor shock associated to our scenario (2)

%data_deep/Data_US_BEA405_DEC2020_demand.xlsx
	This files contains the demand shocks. Here there is none. This is the file where need to paste the user-supplied demand shocks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Matlab Script description

%%%%%%%%%
%%main_BEA504.m
	This is the master file. It run the calibration routine calibrate_BEA504_hat_fun.m, and compute the effect of shocks for various values of elasticities using the routine exercice_hat_fun_nicefigure.m

%%lib/calibrate_BEA504_hat_fun.m
	This is the calibration routine. It takes data from IOUse_Before_Redefinitions_PRO_DET.xlsx (Input-Output tables from the BEA). It computes the value of key parameters. It solve the model for a random labor supply and demand shocks for Cobb-Douglas and user supplied elasticities (or the baseline parameters). To solve the model it uses the system of equation defined in RES_hat_eq.m. This file produce a file in the folder datamatlab/ called US_BEA504_hat.mat. Fianlly it computes network statisitcs and saved them in an excel file.

%%%%%%%%%
%%lib/exercice_hat_fun_nicefigure.m
	This file compute the sectoral effect of shocks. 


	%It takes the following arguments (theta,sigma,epsilon,cap, califile ,excelfile,fig).

		theta,sigma, and epsilon are the elasticities values.

		cap is the maximum value for labor supply shocks.

		califile is the name of the calibration file produced by calibrate_BEA504_hat_fun_bgs.m, its value has to be 'US_BEA504_hat'.

		excelefile has to be equal to the name of the two excel files that contains the laborsupply or demand shocks. Each of this file should have in the first column the list of sector (390 using the BEA classification), and on each column the associated labor (resp. demand) shocks for each sector. The first line should contain a string which will be used for the name of the output file. Each of these files cannot have more than 25 columns.

		fig is equal to 1 if you want to produce and save figures of (log) price/quantity change.


	%Walkthrough 

		It first load the labor_supply shocks from the [excelefile]_laborsupply and the demand shocks from [excelefile]_demand. It combines demand and supply shocks and compute all possible possible pair of demand and labor supply shocks. It also compute scenario with only the labor supply shocks. For each scenario the labor supply and demand shocks are stocked in shocks and demand_shocks respectively.

		It then solves for the system of equation given by RES_hat_eq.m.

		It export in a excel file the response of each aggregate and sectoral variables. The excel file is stored in the folder export. The filename is [califile]_[nameofcolums of supply shocks]_[name of columns of demand shocks if any]_theta[value of theta]_sigma[value of sigma]_epsilon[value of epsilon]_N[value of N].xlsx

		If fig==1 then it also compute a figure of the distribution of (log) quantity across sectors against the data contains in data_deep/Data_quantity.xlsx


%%%%%%%%%
%%lib/RES_hat_eq.m

	This file gives the system of equation that characterized the change of eavery variable following a labor_supply or demand shocks. It need to be given to a solver (e.g. fsolve).





