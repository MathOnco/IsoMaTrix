clc;clear;close all;

%% define colors: (R,G,B);
black = [0,0,0]; red = [1,0,0]; blue = [0,1,0]; green = [0,0,1];

%% label each type (one letter only)
labels = {'A','B','C'};

%% choose payoff matrix:
    % A B C
A =  [0,0,1;  % A
      1,0,0;  % B
      0,1,0]; % C

%% plot total velocity magnitude, with quivers
figure(1); hold on;
title('Total Velocity Magnitude');

% add velocity gradient to the background of simplex:
isomatrix_velocity(A);

% add gridlines to the simplex:
add_gridlines(10); % 10 gridlines in each direction

% add all 3 isoclines, colored in red (1), blue (2) and green (3)
isomatrix_isocline(A);

% add the pairwise equilibriums on the edge of simplex:
isomatrix_fixedpoint(A);

% add the quiver plot:
isomatrix_quiver(A);

% label each corner:
add_labels(labels);



%% velocity of type 1 only
id = 1; % choose type 1, 2, or 3
figure(2); hold on;

% add title
title_str = strcat('Type ', num2str(id),' Velocity');
title(title_str);

% add velocity gradient of only type i's component 
% to the background of simplex:
isomatrix_velocity(A,id);

% add the pairwise equilibriums on the edge of simplex:
isomatrix_fixedpoint(A);

% add isoclines of type i
isomatrix_isocline(A,'LineStyle',':');

% plot some example trajectories (evenly distributed)
isomatrix_trajectory(A);

% label each corner:
add_labels(labels);


%% split domain into regions
figure(3); hold on;
title('IsoMaTrix Regions');

% add the pairwise equilibriums on the edge of simplex:
isomatrix_fixedpoint(A);

% bin state space into regions where 
% sign (+/-) of each type's velocity is
% indicated:
isomatrix_region(A);

% label each corner:
add_labels(labels);
  