function [] = HAL_isomatrix_velocity(id,varargin)

    p = inputParser;
    
    % if odd number of arguments:
    if ((mod(nargin,2) == 0) && (nargin > 1))
        % no user specified id, with extra arguments
        varargin  = [ { id }, varargin ];
        id=[1,2,3];
        addOptional(p,'id',id);
    elseif (nargin == 0)
        % no user specified id, with no extra arguments
        id=[1,2,3];        
    end

    
    %% set up default values for optional parameters: ('Labels')
    labels = {'','',''};
    defaultFilename = 'IsoMaTrixGrid.csv';
        
    % validation of user input labels:
    errorMsg1 = strcat('Labels error: please provide vector of size=3.'); 
    errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
    labelLength = @(x) assert(length(x)==3,errorMsg1);
    labelType = @(x) assert(length(x)==3,errorMsg2); % this is incorrect
    isstringy = @(x) assert(ischar(x),'Filename must be a string');
    addParameter(p,'Labels',labels);
        
    % read in optional parameters
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
            labelLength(labels);
            labelType(labels);
        elseif strcmp(varargin{ind}, 'Filename')
            defaultFilename=varargin{ind+1};
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
        
    %% set up meshgrid
    step = 1/gridlines;
    x_grid = 0:step:(1);
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid); % Generate domain.
    N = sum(x0(1,:));
    grid_step = N/gridlines;
    
    %% remove irrelevant area of domain (outside triangle)
    w = P + Q;
    out = w > 1;
    P(out) = nan;
    Q(out) = nan;
    
    %% coordinate transformation
    y1 = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3); 
    y2 = P * (1/2)*tan(60*pi/180);
    
    %% prepare matrices    
    [rows,~] = size(data);
    
    [n,m] = size(P);
    Z = zeros(n,m);
    Z1 = zeros(n,m);
    Z2 = zeros(n,m);
    Z3 = zeros(n,m);    
    
    for row = 1:sims:rows
        j = round(x0(row,1)/grid_step) + 1;
        i = round(x0(row,2)/grid_step) + 1;
        R = xF(row:(row+sims-1),:) - x0(row:(row+sims-1),:);
        
        Z1(i,j) = mean(R(:,1)); %Z1sd(i,j) = std(R(:,1));
        Z2(i,j) = mean(R(:,2)); %Z2sd(i,j) = std(R(:,2));
        Z3(i,j) = mean(R(:,3)); %Z3sd(i,j) = std(R(:,3));
        
        % mean of the mags, not mag of the means
        Mag = sqrt(R(:,1).^2 + R( :,2).^2 + R( :,3).^2);
        Z(i,j) = mean(Mag); 
        %Zsd(i,j) = std(Mag);
        
        % arrow mag:
%         V = R(:,1); U = (R(:,3)-R(:,2))*cos(pi/3);
%         angle = atan2(V,U);
%         Asd(i,j) = std(angle);
        
    end
    

    if (length(id) == 3)
        BINS = 10;
        
        %% plot total velocity magnitude
        [~,h]=contourf(y1,y2,Z,BINS);hold on;
        set(h,'LineColor','none');
        colormap(gca,parula);
        colorbar('FontSize',16);

    else
        BINS = 40;
        rainbow = [];
        Z = Z3;
%         Zsd = Z3sd;
        if (id == 1)
            Z = Z1;
%             Zsd = Z1sd;
        elseif (id == 2)
            Z = Z2;
%             Zsd = Z2sd;
        end

        mn = min(min(Z));
        mx = max(max(Z));

        blue = [0,0,1];
        white = [1,1,1];
        red = [1,0,0];

        BINS = round(BINS/2);
        
        %% blue to white gradient of colors
        blue_to_white = zeros(BINS-1,3);
        blue_to_white(1,:) = blue;
        blue_to_white(BINS-1,:) = white;
        
        %% white to red gradient of colors
        white_to_red = zeros(BINS-1,3);
        white_to_red(1,:) = white;
        white_to_red(BINS-1,:) = red;

        for i = 2:1:(BINS-2)
            p = (i-1)/(BINS-2);
            
            % blue to white:
            c1 = blue(1) + (white(1) - blue(1)) *p;
            c2 = blue(2) + (white(2) - blue(2)) *p;
            c3 = blue(3) + (white(3) - blue(3)) *p;
            blue_to_white(i,:) = [c1,c2,c3];
            
            % white to red:
            c1 = white(1) + (red(1) - white(1)) *p;
            c2 = white(2) + (red(2) - white(2)) *p;
            c3 = white(3) + (red(3) - white(3)) *p;
            white_to_red(i,:) = [c1,c2,c3];
            
        end
        
        %% this algorithm only works if bluebins = redbins
        rainbow = [blue_to_white(1:end-1,:);
                    white;white;white;
                    white_to_red(2:end,:)];         

        
        my_scale = max(abs(mx),abs(mn));
        
        %% plot uncertainty, too
        colormap(gca,rainbow);
        caxis([-my_scale,my_scale]);
        [~,h]=contourf(y1,y2,Z,BINS);hold on;
        set(h,'LineColor','none');
        colorbar('FontSize',16);

    end
    add_labels(labels);
end
