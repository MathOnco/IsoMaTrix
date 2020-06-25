function [] = isomatrix_region(A,varargin)

    p = inputParser;
    
    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [1,1,1];
    labels = {'','',''};
    linestyle = '-';
    linewidth = 3;
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

    %% velocities
    Z1 = P.*(f1-phi);
    Z2 = Q.*(f2-phi);
    Z3 = (1-P-Q).*(f3-phi);
    Z = Z1;

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


    colormap(gca,rainbow);
    colorbar;
    caxis([1.5,7.5]);
    
    cmapdef = colormap(gca,rainbow) ; %Define Colormap
    cmap = cmapdef(1:8:end, :) ; %Find Values of colors corresponding to each point plotted
    cbh = colorbar('YTickLabel', {'-, +, -','+, +, -','+, -, -','+, -, +','-, -, +','-, +, +'}) ;
    
    cbh.FontSize = 24;
    [~,h]=contourf(y1,y2,Z,100); hold on;
    set(h,'LineColor','none');

    % set isoclines on afterwards
    isomatrix_isocline(A,1,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth);
    isomatrix_isocline(A,2,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth);
    isomatrix_isocline(A,3,'Color',color,'LineStyle',linestyle,'LineWidth',linewidth);
    
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

