function [] = HAL_isomatrix_trajectory(varargin)

    % defaults
    color = [0,0,0];
    labels = {'','',''};
    linestyle = '-';
    linewidth = 2;
    defaultFilename = 'HAL_trajectory.csv';

    p = inputParser;


    
    %% set up default values for optional parameters: ('Color,' 'Labels,' 'LineWidth,' and 'LineStyle')    
    % validation of user input color:
    vectorValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,3]});
    addParameter(p,'Color',color,vectorValidator)
        
    % validation of user input labels:
    errorMsg1 = strcat('Labels error: please provide vector of size=3.'); 
    errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
    labelLength = @(x) assert(length(x)==3,errorMsg1);
    labelType = @(x) assert(length(x)==3,errorMsg2);
    addParameter(p,'Labels',labels);
        
    % validate linestyle, linewidth:
    addParameter(p,'LineStyle',linestyle);
    lwValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,1]});
    addParameter(p,'LineWidth',linewidth,lwValidator);
    isstringy = @(x) assert(ischar(x),'Filename must be a string');
        
    % read in optional parameters    
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Color')
            color=varargin{ind+1};
        elseif strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
            labelLength(labels);
            labelType(labels);
        elseif strcmp(varargin{ind}, 'LineStyle')
            linestyle=varargin{ind+1};
            assert(ischar(linestyle),'Incorrectly specified LineStyle.');
        elseif strcmp(varargin{ind}, 'LineWidth')
            linewidth=varargin{ind+1};
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
    
    %% read in data
    data = dlmread(filepath,',',1,1);
    x_over_time = data(:,1:3);
    x_over_time = x_over_time./sum(x_over_time(1,:));

    %% plot trajector(y/ies)    
    [x_points,y_points] = UVW_to_XY(x_over_time(:,1:3));
    plot(x_points,y_points,linestyle, 'LineWidth', linewidth,'Color',color);hold on;        
    add_labels(labels);

end