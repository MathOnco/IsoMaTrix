function [] = isomatrix_fixedpoint(A,index,varargin)

    p = inputParser;
    % if odd number of arguments:
    if ((mod(nargin,2) == 1) && (nargin > 2))
        % no user specified index, with extra arguments
        varargin  = [ { index }, varargin ];
        index = 1;
        addOptional(p,'index',index);
    elseif (nargin == 1)
        % no user specified index, with no extra arguments
        index = 1;
    end
    
    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [0,0,0];
    labels = {'','',''};
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
        end
    end
               
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    %% fix margins if there are many games
    H = (1/2)*tan(60*pi/180); L = 1;
    del = (1 - H)/2;
    mul = index;
    xlim([-mul*del,L + mul*del]);
    ylim([-mul*del,H+mul*del]);

    % offset axes parameters (probably best not to change)
    line = [0 1];
    REL_DIST = 0.05;
    delta = [0,0];
    
    us_ms = 12;
    s_ms = 50;

    ylim([-REL_DIST*index,(1+REL_DIST*index)/sin(pi/3)]);
    xlim([-REL_DIST*index,1+REL_DIST*index]);

    for i = [1,2]
        for d = [1,2]
            j = mod(i-1 + d,3)+1;
            if (i < j)

                %% calculate offset of two-by-two subgame lines:
                if (i == 2) && (j == 3)
                    delta = [0,-REL_DIST*index];
                elseif (i == 1) && (j == 2)
                    delta = [-REL_DIST*index,REL_DIST*index*tan(pi/6)];
                else
                    delta = [REL_DIST*index,REL_DIST*index*tan(pi/6)]; 
                end

                %% line data for offset line for two-by-two sub-game:
                line_x = zeros(length(line),3);
                line_x(:,i) = line';
                line_x(:,j) = 1-line';
                [x_line,y_line] = UVW_to_XY(line_x);

                %% this is the subgame:
                Ap = [A(i,i), A(i,j);
                         A(j,i), A(j,j)];

                % equal competition
                if ((Ap(1,1) == Ap(2,1)) && (Ap(1,2) == Ap(2,2)))
                    % equal competition game 
                    % ( no interior fixed points)

                    % plot dashed line
                    if (index > 0)
                        plot(x_line + delta(1),y_line + delta(2),':', 'LineWidth', 3,'Color',color);hold on; 
                    end

                % interior point
                elseif ((Ap(1,1)-Ap(2,1))*(Ap(1,2)-Ap(2,2)) < 0)

                    % plot line
                    if (index > 0)
                        plot(x_line + delta(1),y_line + delta(2),'-', 'LineWidth', 3,'Color',color);hold on; 
                    end

                    % there exists an interior point
                    x_star = (Ap(2,2) - Ap(1,2))/ (Ap(1,1) - Ap(1,2) - Ap(2,1) + Ap(2,2) );

                    % [u,v,w] (point is the same if either un/stable)
                    x = zeros(1,3);
                    x(i) = x_star;
                    x(j) = 1 - x_star;
                    [x_point,y_point] = UVW_to_XY(x);

                    if ((Ap(1,1) < Ap(2,1)) && (Ap(1,2) > Ap(2,2)))
                        % this interior point is stable:
                        plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                        plot_arrows(i,j,x_star,delta,color,1);
                        
                    else
                        % this interior point is unstable:
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'Color',[1,1,1],'MarkerEdgeColor',color,'MarkerFaceColor',[1,1,1]);
                        plot_arrows(i,j,x_star,delta,color,-1);
                    end

                else
                    % plot line
                    if (index > 0)
                        plot(x_line + delta(1),y_line + delta(2),'-', 'LineWidth', 3,'Color',color);hold on; 
                    end
                    
                    % testing
                    if (Ap(1,1) >= Ap(2,1))
                        plot_arrows(i,j,1,delta,color,1);
                    else
                        plot_arrows(i,j,1,delta,color,-1);
                    end
                    
                end

                % if not equal, determine other stable point:
                if ((Ap(1,1) == Ap(2,1)) && (Ap(1,2) == Ap(2,2)))
                    
                else
                
                    % determine left stability:
                    x = zeros(1,3);            
                    x(i) = 1;
                    [x_point,y_point] = UVW_to_XY(x);


                    if (Ap(1,1) >= Ap(2,1))
                        plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                    else
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                    end

                    x = zeros(1,3);
                    x(j) = 1;
                    [x_point,y_point] = UVW_to_XY(x);

                    if (Ap(2,2) > Ap(1,2))                    
                        plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                    else
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                    end
                end
                

                %% plot a line over all three edges, just for cleanliness:
                plot(x_line,y_line,'-', 'LineWidth', 3,'Color',[0,0,0]);hold on;


            end

        end  
    end

    %% solve Ax = b for internal equil
    Apq = zeros(2,2);
    Apq(1,1) = A(1,1)-A(1,3)-A(3,1)+A(3,3);
    Apq(1,2) = A(1,2) - A(1,3) - A(3,2) + A(3,3);

    Apq(2,1) = A(2,1)-A(2,3)-A(3,1)+A(3,3);
    Apq(2,2) = A(2,2)-A(2,3)-A(3,2)+A(3,3);

    B = zeros(2,1);
    B(1) = A(3,3) - A(1,3);
    B(2) = A(3,3) - A(2,3);

    pq_star = inv(Apq)*B;
    p_star = pq_star(1);
    q_star = pq_star(2);
    
    if ((p_star > 0) && (q_star > 0))
        if ((p_star < 1) && (q_star < 1))
            if (((p_star+q_star) < 1))
                x = [p_star,q_star,1-p_star-q_star];
                
                type = DetermineFixedPointType(x,A);
                [x_point,y_point] = UVW_to_XY(x);
                if ((type <= 1))
                    % stable:
                    plot(x_point,y_point,'.', 'MarkerSize', s_ms,'Color',color);hold on;
                elseif (type == 2)
                    % unstable
                    plot(x_point,y_point,'^','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);hold on;
                elseif (type == 3)
                    % saddle (triangle)
                    plot(x_point,y_point,'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);hold on;
                else
                    % source (open circle)
                    plot(x_point,y_point,'s','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);hold on;
                end
            end
        end
    end
    
    add_labels(labels);
end

function plot_arrows(i,j,x_star,delta,color,stability)
    % direction:
    xDiff = zeros(1,3);
    xDiff(i) = 1;
    xDiff(j) = -1;

    % down-left
    V = stability*xDiff(1);
    U = stability*(xDiff(3)-xDiff(2))* cos(pi/3);

    arrow_length = -0.03;
    
    if (stability < 0)
        arrow_length = 0.03;
    end
    
    
    len = 0.00001;
    
    xMid = zeros(1,3);
    xMid(i) = x_star/2-arrow_length;
    xMid(j) = 1 - x_star/2+arrow_length;
    [xM,yM] = UVW_to_XY(xMid);
    
    hw = 16;
    hl = 20;
    
    %% bottom up (ignore if too close to edge)
    if (x_star > 0.1)
        ah = annotation('arrow','headStyle','cback1','HeadLength',hl,'HeadWidth',hw);
        set(ah,'parent',gca);
        set(ah,'position',[xM+delta(1),yM+delta(2), len*U, len*V]);
        set(ah,'Color',color);
    end
    
    %% top down (ignore if too close to edge)
    if (x_star < 0.9)
        xMid = zeros(1,3);
        xMid(i) = 1-(1-x_star)/2+arrow_length;
        xMid(j) = (1-x_star)/2-arrow_length;
                
        [xM,yM] = UVW_to_XY(xMid);
        
        ah = annotation('arrow','headStyle','cback1','HeadLength',hl,'HeadWidth',hw);
        set(ah,'parent',gca);
        set(ah,'position',[xM+delta(1),yM+delta(2), -len*U, -len*V]);
        set(ah,'Color',color);        
    end
end


function [type] = DetermineFixedPointType(x,A)

    [~,lambda,~]=hessian(x,A);

    lambda1 = real(lambda(1,1));
    lambda2 = real(lambda(2,2));
    
    if ((lambda1== 0) && (lambda2==0))
        % this is stable, but not asymptotically stable
        type = 0;
    elseif (((lambda1== 0) && (lambda2>0)) || ((lambda1>0) && (lambda2==0)) )
        % one zero, one negative
        type = 0;
    elseif (((lambda1== 0) && (lambda2<0)) || ((lambda1< 0) && (lambda2==0)) )
        % one zero, one positive == UNSTABLE
        type = 2;
    elseif lambda1*lambda2<0   
        % this is a saddle point
        type = 3;
    elseif ((lambda1> 0) && (lambda2>0))
        % this is a source
        type = 4;
    else
        % this is an attractor
        type = 1;
    end
end


