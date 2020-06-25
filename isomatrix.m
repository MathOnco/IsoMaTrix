function [] = isomatrix(A,varargin)


p = inputParser; 
labels = {'1','2','3'};
[nn,~] = size(A);

assert(3==nn,'Please provide a 3 by 3 matrix.')

% validation of user input labels:
errorMsg1 = strcat('Labels error: please provide vector of size=',' ',num2str(nn),').'); 
errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
labelLength = @(x) assert(length(x)==nn,errorMsg1);
labelType = @(x) assert(length(x)==nn,errorMsg2);
addParameter(p,'Labels',labels);

% read in optional parameters
[nParams] = length(varargin);
for param = 1:1:(nParams/2)
    index = (param-1)*2 + 1;
    if strcmp(varargin{index}, 'Labels')
        labels=varargin{index+1};
        labelLength(labels);
        labelType(labels);
    end
end


%% plot total velocity magnitude, with quivers
figure(1); hold on;
title('Total Velocity Magnitude');
isomatrix_velocity(A,1:3);
isomatrix_fixedpoint(A,'Color',[0,0,0]);
isomatrix_quiver(A,'Color',[0,0,0]);
add_labels(labels);

%% velocity of type 1 only

for id = 1:3
    figure(id + 1); hold on;
    title_str = strcat('Type-', num2str(id),' Velocity');
    title(title_str);
    isomatrix_velocity(A,id);
    isomatrix_fixedpoint(A,'Color',[0,0,0]);
    isomatrix_isocline(A,id,'Color',[0,0,0],'LineStyle',':');
    add_labels(labels);
end

%% split domain into regions
figure(5); hold on;
title('IsoMaTrix Regions');
isomatrix_fixedpoint(A,'Color',[0,0,0]);
isomatrix_region(A);
add_labels(labels);

end