function [] = HAL_isomatrix_region(varargin)

    labels = {'','',''};
    defaultFilename = 'IsoMaTrixGrid.csv';
    
    p = inputParser;
    addParameter(p,'Labels',labels);
    
    % validation of user input labels:
    errorMsg1 = strcat('Labels error: please provide vector of size=3.'); 
    errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
    isstringy = @(x) assert(ischar(x),'Filename must be a string');
    labelLength = @(x) assert(length(x)==3,errorMsg1);
    labelType = @(x) assert(length(x)==3,errorMsg2); % this is incorrect

    % read in optional parameters
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        index = (param-1)*2 + 1;
        if strcmp(varargin{index}, 'Labels')
            labels=varargin{index+1};
            labelLength(labels);
            labelType(labels);
        elseif strcmp(varargin{index}, 'Filename')
            defaultFilename=varargin{index+1};
            % check if string
            isstringy(defaultFilename);
        end
    end

    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    % string concatenate filepath:
    filepath = strcat('HALMatrixGame/HALMatrix-output/',defaultFilename);

    %% read in data
    data = dlmread(filepath,',',2,0);
    x0 = data(:,1:3);
    xF = data(:,4:6);
    sims = max(data(:,7))+1;
    gridlines = max(data(:,8));
        
    %% create mesh grid
    step = 1/gridlines;
    x_grid = 0:step:(1);
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid); % Generate domain.
    
    %% remove irrelevant:
    w = P + Q;
    out = w > 1;
    P(out) = nan;
    Q(out) = nan;
    
    % coordinate transformation
    y1 = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3); 
    y2 = P * (1/2)*tan(60*pi/180);
    
    [rows,~] = size(data);
    
    [n,m] = size(P);
    Z = zeros(n,m);
    Z1 = zeros(n,m);
    Z2 = zeros(n,m);
    Z3 = zeros(n,m);
        
    N = sum(x0(1,:));
    
    grid_step = N/gridlines;
    
    for row = 1:sims:rows
        j = round(x0(row,1)/grid_step) + 1;
        i = round(x0(row,2)/grid_step) + 1;
        R = xF(row:(row+sims-1),:) - x0(row:(row+sims-1),:);
        
        Z1(i,j) = mean(R(:,1));
        Z2(i,j) = mean(R(:,2));
        Z3(i,j) = mean(R(:,3));
        
        % mean of the mags, not mag of the means
        Mag = sqrt(R(:,1).^2 + R( :,2).^2 + R( :,3).^2);
        Z(i,j) = mean(Mag); 
        
        % arrow mag:
        V = R(:,1); U = (R(:,3)-R(:,2))*cos(pi/3);
        angle = atan2(V,U);
        
    end

    blue = [0.2188,0.4531,0.6914];
    red = [0.7734,0.2188,0.1719];
    green = [0.3086,0.6211,0.2227];
    blue_red = [229,148,159]/255;
    blue_green = [170,229,235]/255;
    red_green = [208,184,74]/255;

    rainbow = [green;
        red_green;
        red;
        blue_red;
        blue;
        blue_green];
    
    [n,m] = size(Z);

    for i = 1:n
        for j = 1:m
            if (~isnan(Z(i,j)))
                %% test all cases
                if ( (Z1(i,j)>=0) && (Z2(i,j)>=0) && (Z3(i,j)>=0) )
                    Z(i,j) = 1; % 000
                elseif ( (Z1(i,j)>=0) && (Z2(i,j)>=0) && (Z3(i,j)<0) )
                    Z(i,j) = 3; % 001
                elseif ( (Z1(i,j)>=0) && (Z2(i,j)<0) && (Z3(i,j)>=0) )
                    Z(i,j) = 5; % 010
                elseif ( (Z1(i,j)>=0) && (Z2(i,j)<0) && (Z3(i,j)<0) )
                    Z(i,j) = 4; % 011
                elseif ( (Z1(i,j)<0) && (Z2(i,j)>=0) && (Z3(i,j)>=0) )
                    Z(i,j) = 7; % 100
                elseif ( (Z1(i,j)<0) && (Z2(i,j)>=0) && (Z3(i,j)<0) )
                    Z(i,j) = 2; % 101
                elseif ( (Z1(i,j)<0) && (Z2(i,j)<0) && (Z3(i,j)>=0) )
                    Z(i,j) = 6; % 110
                elseif ( (Z1(i,j)<0) && (Z2(i,j)<0) && (Z3(i,j)<0) )
                    Z(i,j) = 8; % 111
                else
                    Z(i,j) = -1;
                end

            end
        end
    end   

    gcf; hold on;
    colormap(gca,rainbow);
    colorbar;
    caxis([1.5,7.5]);
    
    cmapdef = colormap(gca,rainbow) ; %Define Colormap
    cmap = cmapdef(1:8:end, :) ; %Find Values of colors corresponding to each point plotted
    cbh = colorbar('YTickLabel', {'-, +, -','+, +, -','+, -, -','+, -, +','-, -, +','-, +, +'}) ;
    cbh.FontSize = 24;
    
    
    [~,h]=contourf(y1,y2,Z,100); hold on;
    set(h,'LineColor','none');
    
    L = 1;
    H = (L/2)*tan(60*pi/180);
    d1 = 0.26;
    d2 = 0.18;
    
    legend_string = {'($\dot{x}_1, \dot{x}_2, \dot{x}_3$)'};
    if ~isempty(labels{1}) 
         legend_string = {strcat('($\dot{x}_',labels{1}, ', \dot{x}_',labels{2}, ', \dot{x}_',labels{3},'$)')};
    end
    
    text(L+d1,H+d2,legend_string, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 18, 'Interpreter','Latex');

    

    add_labels(labels);


end

