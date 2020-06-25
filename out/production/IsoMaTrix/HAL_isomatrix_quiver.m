function [] = HAL_isomatrix_quiver(uncertainty_boolean,varargin)
        

    p = inputParser;
    
    % if odd number of arguments:
    if ((mod(nargin,2) == 0) && (nargin > 1))
        % no user specified id, with extra arguments
        varargin  = [ { uncertainty_boolean }, varargin ];
        uncertainty_boolean = true;
        addOptional(p,'uncertainty_boolean',uncertainty_boolean);
        
    elseif (nargin == 0)
        % no user specified uncertainty (default is on)
        uncertainty_boolean = true;       
    end

    
    %% set up default values for optional parameters: ('Labels', 'Color','Filename')
    labels = {'','',''};
    defaultFilename = 'IsoMaTrixGrid.csv';
    color = [0,0,0];
    
    vectorValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,3]});
    addParameter(p,'Color',color,vectorValidator)
        
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
        elseif strcmp(varargin{ind}, 'Color')
            color=varargin{ind+1};
        elseif strcmp(varargin{ind}, 'Filename')
            defaultFilename=varargin{ind+1};
            isstringy(defaultFilename); % check if string
        end
    end
    
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;

    
    % string concatenate filepath:
    filepath = strcat('HALMatrixGame/HALMatrix-output/',defaultFilename);


    if (~exist('color','var') || isempty(color))
        color = [0,0,0]; % black arrows if none specified
    end

    if (~exist('filepath','var') || isempty(filepath))
        filepath  = 'HALMatrixGame/HALMatrix-output/IsoMaTrixGrid.csv';
    end
    
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    %% read in data from HALMatrix
    data = dlmread(filepath,',',2,0);
    x0 = data(:,1:3);
    xF = data(:,4:6);
    sims = max(data(:,7))+1;
    gridlines = max(data(:,8));
    N = sum(x0(1,:));
    grid_step = N/gridlines;
       
    %% set up mesh grid
    step = 1/gridlines;
    x_grid = 0:step:(1);
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid); % Generate domain.
    
    %% remove irrelevant:
    w = P + Q;
    out = w > 1;
    P(out) = nan;
    Q(out) = nan;
        
    [rows,~] = size(data);
    
    [n,m] = size(P);
    Z = zeros(n,m); Zsd = zeros(n,m);
    Z1 = zeros(n,m); Z1sd = zeros(n,m);
    Z2 = zeros(n,m); Z2sd = zeros(n,m);
    Z3 = zeros(n,m); Z3sd = zeros(n,m);
    
    % uncertainty of arrow direction:
    Asd = zeros(n,m);  
    
    
   
    
    for row = 1:sims:rows
        j = round(x0(row,1)/grid_step) + 1;
        i = round(x0(row,2)/grid_step) + 1;
        R = xF(row:(row+sims-1),:) - x0(row:(row+sims-1),:);
        
        
        Z1(i,j) = mean(R(:,1)); 
        Z1sd(i,j) = std(R(:,1));
        Z2(i,j) = mean(R(:,2)); 
        Z2sd(i,j) = std(R(:,2));
        Z3(i,j) = mean(R(:,3)); 
        Z3sd(i,j) = std(R(:,3));
        
        % mean of the mags, not mag of the means
        Mag = sqrt(R(:,1).^2 + R( :,2).^2 + R(:,3).^2);
        Z(i,j) = mean(Mag); 
        Zsd(i,j) = std(Mag);
        
        % arrow mag:
        V = R(:,1); 
        U = (R(:,3)-R(:,2))*cos(pi/3);
        angle = atan2(V,U);
        Asd(i,j) = std(angle);
        
    end
    
    

    %% arrow length:
    len = 0.0375/gridlines*20;
    HL = 8;
    HW = 8;
    
    step0 =step/2;

    %% coordinate transformation:
    X = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3);
    Y = P * (1/2)*tan(60*pi/180);
        

    %% quiver
    for i = 1:length(P)
        for j = 1:length(Q)                
            if ( P(i,j) < step0) || (Q(i,j) < step0) || ( (P(i,j)+Q(i,j)) > 1-step0)
                %% outside of simplex

            else

                %% calculate direction & magnitude
                V = Z1(i,j);
                U = (Z3(i,j)-Z2(i,j))*cos(pi/3); 
                mag = sqrt(V^2 + U^2);

                if (mag > 0)


                    %% plot uncertainty arc:
                    if uncertainty_boolean
                        angle = atan2(len*V/mag,len*U/mag);
                        a = angle - Asd(i,j);
                        b = angle + Asd(i,j);
                        h=UncertaintyArc(a,b,X(i,j), Y(i,j),len*0.5,color);hold on;
                        set(h,'edgecolor','none','facecolor',color,'facealpha',0.5);
                    end

                    %% set direction & magnitude was 8 and 5
                    ah = annotation('arrow','headStyle','cback1','HeadLength',HL,'HeadWidth',HW);
                    set(ah,'parent',gca);
                    set(ah,'position',[X(i,j), Y(i,j), len*U/mag, len*V/mag]);
                    set(ah,'Color',color);


                end
            end
        end
    end
    add_labels({'','',''});
end


function P = UncertaintyArc(a1,a2,height,offset,radius,color)
x = radius*cos(linspace(a1,a2)) + height;
y = radius*sin(linspace(a1,a2)) + offset;
x = [x height x(1)];
y = [y offset y(1)];
P = fill(x,y,color);
end
