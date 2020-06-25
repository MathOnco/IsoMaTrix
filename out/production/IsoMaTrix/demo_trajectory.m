clc;clear;close all;

%% define colors: (R,G,B);
black = [0,0,0]; red = [1,0,0];

%% label each type (one letter only)
labels = {'A','B','C'};

%% choose payoff matrix:
    % A B C
A =  [0,0,1;  % A
      1,0,0;  % B
      0,1,0]; % C
 
%% plot total velocity magnitude, with quivers
figure(1); hold on;
isomatrix_velocity(A);
colorbar;

add_gridlines(10); % 10 gridlines in each direction
isomatrix_fixedpoint(A);
isomatrix_quiver(A);
add_labels(labels);

%% plot a trajectory on simplex starting at x0 until time=tF
x0 = [0.1,0.2,0.7];
tF = 10; 
isomatrix_trajectory(A,x0,tF,'Color',red);

%% line plot of same simulation
figure(2);
line_plot(A,x0,tF,'Labels',labels);
