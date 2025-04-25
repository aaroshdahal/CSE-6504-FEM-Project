4%%                           MATLAB PROJECT                  
%__________________________________________________________________________
% 
%                    CEE-6504 Finite Element Methods
%                            Dr. Rafi Muhanna
%                   Georgia Institute of Technology
%                              Spring 2025                     
%__________________________________________________________________________
% 
%                    Developed by  Shahrokh Shahi 
%__________________________________________________________________________
%
%
%% Initialization (Do NOT modify this section)

clc
clear
close all

format short g
format compact

% adding "lib" folders to MATLAB path
path = mfilename('fullpath');
path(end-length(mfilename):end)=[];
addpath(fullfile(path,'lib'));
addpath(fullfile(path,'lib_io'));

%% input file name

% you can input the data file name here:
inpFileName = 'input2_beamlike.txt';

%% Reading the Input File

% interprete the input file and storing all data in the "Model" structure
Model = inpFileReader(inpFileName);

% displaying a summary of the input data (Model structure) on the Command
% Window
 printSummary(Model);

%% Initial Visualization

plotMeshX(Model)

%% Call the FEM function

% TODO (just uncomment this line, whenever the fem2D.m becomes completed)
Results = fem2D(Model);
printOutput(Results)
