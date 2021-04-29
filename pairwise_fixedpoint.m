function [] = pairwise_fixedpoint(A,varargin)

    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [0,0,0];
    labels = {};
    [nn,~] = size(A);
    for l = 1:nn
        labels{l} = num2str(l);
    end
    
    p = inputParser;
    
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
        index = (param-1)*2 + 1;
        if strcmp(varargin{index}, 'Color')
            color=varargin{index+1};
        elseif strcmp(varargin{index}, 'Labels')
            labels=varargin{index+1};
            labelLength(labels);
            labelType(labels);
        end
    end
           
    h = gcf;clf;
    figure_number=h.Number;
    figure(figure_number); hold on;

    % offset axes parameters
    line = [0 1];
    us_ms = 12;
    s_ms = 50;
    REL_DIST = 0.05;
    bump_index = 1;
    
    ivec = 1:1:(nn-1);
    
    for i = ivec
        for d = ivec
            j = mod(i-1 + d,nn)+1;
            if (i < j)
                                
                ii = 2;
                jj = 3;

                delta = [0,-REL_DIST*bump_index];
                line_x = zeros(length(line),3);
                line_x(:,ii) = line';
                line_x(:,jj) = 1-line';
                [x_line,y_line] = UVW_to_XY(line_x);
                
                %% this is the subgame:
                Ap = [A(i,i), A(i,j); A(j,i), A(j,j)];
                
                %% add gradient of selection:
                xfit=0:0.001:1;
                [~,xdot2,~] = find_xdot(xfit,Ap);
                MAX_MIN = max(abs(min(xdot2)),max(xdot2));                
                myNorm = REL_DIST*1000*MAX_MIN;
                plot(xfit + delta(1),xdot2./myNorm + delta(2),'-', 'LineWidth', 1.5,'Color',color);hold on;
                
                % color area beneath the gradient curve:
                xx = xfit + delta(1);                
                curve1 = xx.*0 + delta(2);
                curve2 = xdot2./myNorm + delta(2);
                inBetween = [curve1, fliplr(curve2)];
                x2 = [xx, fliplr(xx)];
                fill(x2, inBetween, color/2,'FaceAlpha',0.2,'LineStyle','none');
                
                % label left / right
                del = 0.04;
                text(-del,y_line(1)+ delta(2),labels{i}, 'HorizontalAlignment','right', 'VerticalAlignment','middle', 'fontsize', 18 );
                text(1+del,y_line(1)+ delta(2),labels{j}, 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 18 );

                

                % equal competition
                if ((Ap(1,1) == Ap(2,1)) && (Ap(1,2) == Ap(2,2)))
                    % equal competition game 
                    % ( no interior fixed points)

                    % plot dashed line
                    plot(x_line + delta(1),y_line + delta(2),':', 'LineWidth', 3,'Color',color);hold on; 

                % interior point
                elseif ((Ap(1,1)-Ap(2,1))*(Ap(1,2)-Ap(2,2)) < 0)

                    % plot line
                    plot(x_line + delta(1),y_line + delta(2),'-', 'LineWidth', 3,'Color',color);hold on; 

                    % there exists an interior point
                    x_star = (Ap(2,2) - Ap(1,2))/ (Ap(1,1) - Ap(1,2) - Ap(2,1) + Ap(2,2) );

                    % [u,v,w] (point is the same if either un/stable)
                    x = zeros(1,3);
                    x(ii) = x_star;
                    x(jj) = 1 - x_star;
                    [x_point,y_point] = UVW_to_XY(x);

                    if ((Ap(1,1) < Ap(2,1)) && (Ap(1,2) > Ap(2,2)))
                        % this interior point is stable:
                        plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                        plot_arrows(ii,jj,x_star,delta,color,1);
                        
                    else
                        % this interior point is unstable:
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'Color',[1,1,1],'MarkerEdgeColor',color,'MarkerFaceColor',[1,1,1]);
                        plot_arrows(ii,jj,x_star,delta,color,-1);
                    end

                else
                    % plot line
                    plot(x_line + delta(1),y_line + delta(2),'-', 'LineWidth', 3,'Color',color);hold on; 
                    
                    % testing
                    if (Ap(1,1) >= Ap(2,1))
                        plot_arrows(ii,jj,1,delta,color,1);
                    else
                        plot_arrows(ii,jj,1,delta,color,-1);
                    end
                    
                end

                % determine left stability:
                x = zeros(1,3);            
                x(ii) = 1;
                [x_point,y_point] = UVW_to_XY(x);


                if (Ap(1,1) >= Ap(2,1))
                    plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                else
                    plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                end

                x = zeros(1,3);
                x(jj) = 1;
                [x_point,y_point] = UVW_to_XY(x);

                if (Ap(2,2) >= Ap(1,2))
                    plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                else
                    plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                end

                
                bump_index = bump_index+1;
            end
        end  
    end
    
    %% clean up domain:
    xlim([-REL_DIST,1+REL_DIST]);
    ylim([-bump_index*REL_DIST,0]);
    set(gca,'visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on'); % leave on title
    box off;
    set(gcf,'color','w');
    width = 600;
    set(gcf,'Position',[100,100,width/sin(pi/3),width]);

    
end

function plot_arrows(i,j,x_star,delta,color,stability)
    % direction:
    xDiff = zeros(1,3);
    xDiff(i) = 1;
    xDiff(j) = -1;

    % down-left

    V = stability*xDiff(1);
    U = stability*(xDiff(3)-xDiff(2))* cos(pi/3);

    arrow_length = 0;
    len = 0.00001;
    
    xMid = zeros(1,3);
    xMid(i) = x_star/2-arrow_length;
    xMid(j) = 1 - x_star/2+arrow_length;
    [xM,yM] = UVW_to_XY(xMid);
    
    hw = 16;
    hl = 20;
    
    %% bottom up
    if (x_star > 0.1)
        ah = annotation('arrow','headStyle','cback1','HeadLength',hl,'HeadWidth',hw);
        set(ah,'parent',gca);
        set(ah,'position',[xM+delta(1),yM+delta(2), len*U, len*V]);
        set(ah,'Color',color);
    end
    
    if (x_star < 0.9)
        xMid = zeros(1,3);
        xMid(i) = 1-(1-x_star)/2+arrow_length;
        xMid(j) = (1-x_star)/2;
                
        [xM,yM] = UVW_to_XY(xMid);
        
        ah = annotation('arrow','headStyle','cback1','HeadLength',hl,'HeadWidth',hw);
        set(ah,'parent',gca);
        set(ah,'position',[xM+delta(1),yM+delta(2), -len*U, -len*V]);
        set(ah,'Color',color);        
    end

end



%% adding this (2 by 2) xdot function:
function [xdot1,xdot2,phi] = find_xdot(x,A)   
    f1 = A(1,1).*(1-x) + A(1,2).*x;
    f2 = A(2,1).*(1-x) + A(2,2).*x;
    phi = f1.*(1-x) + f2.*x;
    
    xdot1 = (1-x).* ( f1 - phi ) ; 
    xdot2 = x.* ( f2 - phi ) ;    
end


