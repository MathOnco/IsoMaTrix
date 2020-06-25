function [] = HAL_isomatrix(varargin)

    % optional arguments: 'Labels' and Fiilename
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
        elseif strcmp(varargin{ind}, 'Filename')
            defaultFilename=varargin{ind+1};
            % check if string
            isstringy(defaultFilename);
        end
    end


    figure(1); hold on;
    title('HALMatrix');
    HAL_isomatrix_velocity('Filename',defaultFilename,'Labels',labels);
    HAL_isomatrix_quiver(false,'Filename',defaultFilename,'Labels',labels);
    add_labels(labels);

    figure(2); hold on;
    title('HALMatrix Uncertainty');
    HAL_isomatrix_uncertainty('Filename',defaultFilename,'Labels',labels);
    HAL_isomatrix_quiver(true,'Filename',defaultFilename,'Labels',labels);
    add_labels(labels);

    %% plot total velocity magnitude, with quivers
    id = 1;
    for f = [3,4,5]
        figure(f); hold on;
        title_str = strcat('Type-', num2str(id),' Velocity');
        title(title_str);
        HAL_isomatrix_velocity(id,'Filename',defaultFilename,'Labels',labels);
        add_labels(labels);
        id = id+1;
    end

    %% split domain into regions
    figure(6); hold on;
    title('IsoMaTrix Regions');
    HAL_isomatrix_region('Filename',defaultFilename,'Labels',labels);
    add_labels(labels);

end