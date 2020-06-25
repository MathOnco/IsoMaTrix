function [] = isomatrix_isocline(A,id,varargin)

    RBG = [1,0,0;0,1,0;0,0,1];

    p = inputParser;
        
    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [1,0,0];
    labels = {'','',''};
    linestyle = '-';
    linewidth = 2;
    [nn,~] = size(A);
    assert(3==nn,'Please provide a 3 by 3 matrix.')
    
    % validation of user input color:
    vectorValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,3]});
    addParameter(p,'Color',color,vectorValidator)
        
    % validation of user input labels:
    errorMsg1 = strcat('Labels error: please provide vector of size=',' ',num2str(nn),').'); 
    errorMsg2 = 'Incorrect label formatting (must be cell-array).'; 
    labelLength = @(x) assert(length(x)==nn,errorMsg1);
    labelType = @(x) assert(length(x)==nn,errorMsg2);
    addParameter(p,'Labels',labels);
    
    % validate linestyle, linewidth:
    addParameter(p,'LineStyle',linestyle);
    lwValidator = @(x) validateattributes(x,{'numeric'},{'size',[1,1]});
    addParameter(p,'LineWidth',linewidth,lwValidator);
    
    if ((mod(nargin,2) == 1) && (nargin > 2))
        varargin  = [ { id },  varargin ];        
        % no user specified id
        id = [1,2,3];
        addOptional(p,'id',id);
    elseif (nargin == 1)
        id = [1,2,3];
    end
    
        
    % read in optional parameters    
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Color')
            color=varargin{ind+1};
            RBG = [color;color;color];
        elseif strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
            labelLength(labels);
            labelType(labels);
        elseif strcmp(varargin{ind}, 'LineStyle')
            linestyle=varargin{ind+1};
            assert(ischar(linestyle),'Incorrectly specified LineStyle.');
        elseif strcmp(varargin{ind}, 'LineWidth')
            linewidth=varargin{ind+1};
        end
    end
               
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;

    %% construct mesh:
    step = 0.005;
    x_grid = 0:step:1;
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid); % Generate domain.

    %% get ride of nullspace in the mesh:
    w = P + Q;
    out = w > 1;
    P(out) = nan;
    Q(out) = nan;

    %% size of triangle (in drawing units):
    y1 = (P./(tan(pi/3)) + (1-P-Q)./(sin(pi/3)))*sin(pi/3); 
    y2 = P * (1/2)*tan(60*pi/180);
    
    %% replicator equation dynamics:
    f1 = (A(1,1) - A(1,3) ).*P + (A(1,2) - A(1,3)).*Q + A(1,3);
    f2 = (A(2,1) - A(2,3) ).*P + (A(2,2) - A(2,3)).*Q + A(2,3);
    f3 = (A(3,1) - A(3,3) ).*P + (A(3,2) - A(3,3)).*Q + A(3,3);
    phi = P.*f1 + Q.*f2 + (1 - P - Q).*f3;
    
    
    
    for strategy = id
        
        Z = (1-P-Q).*(f3-phi);
        if strategy == 1
            Z = P.*(f1-phi);
        elseif strategy == 2
            Z = Q.*(f2-phi);
        end

        if (length(id) == 3)
            contour(y1,y2,Z,[0,0],linestyle,'LineWidth',linewidth,'Color',RBG(strategy,:)); hold on;
        else
            contour(y1,y2,Z,[0,0],linestyle,'LineWidth',linewidth,'Color',color); 
        end
    end
    add_labels(labels);
end