clc;clear;close all;

%% example payoff matrix:

    % A   B   C
A = [0.7,0.0,0.7;  % A
     0.3,0.4,0.8;  % B
     1.0,0.3,0.2]; % C

%% "isomatrix" function creates the following plots:
% - main plot of velocities with quiver velocity gradient
% - velocity component (dx/dt) for each type, with isoclines
% - "isometric _regions" with sign (+/-) for each type's velocity

isomatrix(A);

% or, add labels to each of these plots:
% labels = {'1','2','3'};
% isomatrix(A,labels);









