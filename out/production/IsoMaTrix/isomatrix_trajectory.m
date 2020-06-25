function [] = isomatrix_trajectory(A,x0,tF,varargin)

    % 4 valid calls:
    
    % 1. isomatrix_trajectory(A)
    % 2. isomatrix_trajectory(A,x0,tF)
    % 3. isomatrix_trajectory(A,x0,tF, ... varargins)
    % 4. isomatrix_trajectory(A, ... varargins)
    
    % defaults
    color = [0,0,0];
    labels = {'','',''};
    linestyle = '-';
    linewidth = 2;

    p = inputParser;
    % if odd number of arguments:
    assert(mod(nargin,2)==1,'Please provide correct input arguments: A, x0, tF.');

    if ((nargin > 1) && (ischar(x0)))
        varargin  = [ { x0 }, { tF },  varargin ];        
        % no user specified x0, tF
        tF = 50; % run for 50 timesteps if none specified
        x0 = even_x0();
        addOptional(p,'tF',tF);
        addOptional(p,'x0',x0);
    elseif (nargin == 1)
        % user specified A only
        tF = 50; % run for 50 timesteps if none specified
        x0 = even_x0();
        addOptional(p,'tF',tF);
        addOptional(p,'x0',x0);
    end

    
    %% set up default values for optional parameters: ('Color' and 'Labels')
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
        end
    end
               
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;


    %% plot trajector(y/ies)
    [n,~] = size(x0);
    for i = 1:1:n
        x00 = x0(i,:);        
        [~, xx]=ode45(@(t,n)replicator(t,n,A), [0 tF], x00);
        [x_points,y_points] = UVW_to_XY(xx(:,1:3));
        plot(x_points,y_points,linestyle, 'LineWidth', linewidth,'Color',color);hold on;        
    end

    add_labels({'','',''});

end



function x0 = even_x0()
    step = 0.1;
    step0 =step/2;
    x_grid = 0:step:(1+step); 
    y_grid = x_grid;
    [P,Q] = meshgrid(x_grid,y_grid);
    x0 = zeros(1,3);
    index = 1;
    for i = 1:length(P)
        for j = 1:length(Q)                
            if ( P(i,j) < step0) || (Q(i,j) < step0) || ( (P(i,j)+Q(i,j)) > 1-step0)
                %% outside of simplex
            else
               x0(index,:) = [P(i,j),Q(i,j), 1-P(i,j)-Q(i,j)];
               index = index + 1;
            end
        end
    end
end