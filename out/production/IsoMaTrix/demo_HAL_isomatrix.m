clc;clear;close all;

%% "isomatrix" function creates the following plots, derived from simulations
% of HALMatrixGames. 

%% Data from these simulations are stored previously in:
% 1. HALMatrixGame/HALMatrix-output/HAL_trajectory.csv
% 2. HALMatrixGame/HALMatrix-output/IsoMaTrixGrid.csv

%% this will display
% - main plot of velocities with quiver velocity gradient
% - uncertainty of velocity magnitude
% - velocity component (dx/dt) for each type, with isoclines
% - "isometric _regions" with sign (+/-) for each type's velocity
% - corresponding isomatrix diagram for replicator dynamics w/ identical
% payoff matrix

HAL_isomatrix();







