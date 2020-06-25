function [] = isomatrix_velocity(A,id,varargin)

    p = inputParser;
    % if odd number of arguments:
    if ((mod(nargin,2) == 1) && (nargin > 2))
        % no user specified id, with extra arguments
        varargin  = [ { id }, varargin ];
        id=[1,2,3];
        addOptional(p,'id',id);
    elseif (nargin == 1)
        % no user specified id, with no extra arguments
        id=[1,2,3];
    end

    
    %% set up default values for optional parameters: ('Labels')
    labels = {'','',''};
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
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
            labelLength(labels);
            labelType(labels);
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
    Z1 = P.*(f1-phi);
    Z2 = Q.*(f2-phi);
    Z3 = (1 - P - Q).*(f3-phi);

    if (length(id) == 3)
        BINS = 10;
        %% plot total velocity magnitude
        Z = sqrt(Z1.^2 + Z2.^2 + Z3.^2);
        [~,h]=contourf(y1,y2,Z,BINS);
        set(h,'LineColor','none');
        colormap(gca,parula);
        colorbar('FontSize',16);
    else
        BINS = 40;
        rainbow = [];
        Z = Z3;
        if (id == 1)
            Z = Z1;
        elseif (id == 2)
            Z = Z2;
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
        
        colormap(gca,rainbow);
        caxis([-my_scale,my_scale]);
        [~,h]=contourf(y1,y2,Z,BINS);
        set(h,'LineColor','none');
        colorbar('FontSize',16);
        
    end
    add_labels(labels);
    
end


