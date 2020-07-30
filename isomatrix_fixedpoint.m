function [] = isomatrix_fixedpoint(A,index,varargin)

    p = inputParser;
    effective_index = 1;
    % if odd number of arguments:
    if ((mod(nargin,2) == 1) && (nargin > 2))
        % no user specified index, with extra arguments
        varargin  = [ { index }, varargin ];
        index = 1;
        effective_index = 0;
        addOptional(p,'index',index);
    elseif (nargin == 1)
        % no user specified index, with no extra arguments
        index = 1;
        effective_index = 0;
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
                        
                        CustomMark(x_point,y_point,x,A,color,effective_index);
                        
                    else
                        % this interior point is unstable:
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'Color',[1,1,1],'MarkerEdgeColor',color,'MarkerFaceColor',[1,1,1]);
                        plot_arrows(i,j,x_star,delta,color,-1);
                        
                        CustomMark(x_point,y_point,x,A,color,effective_index);

                    end

                else
                    % plot line
                    if (index > 0)
                        plot(x_line + delta(1),y_line + delta(2),'-', 'LineWidth', 3,'Color',color);hold on; 
                    end
                    
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
                        CustomMark(x_point,y_point,x,A,color,effective_index);
                    else
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                        CustomMark(x_point,y_point,x,A,color,effective_index);
                    end

                    x = zeros(1,3);
                    x(j) = 1;
                    [x_point,y_point] = UVW_to_XY(x);

                    if (Ap(2,2) > Ap(1,2))                    
                        plot(x_point+delta(1),y_point+delta(2),'.', 'MarkerSize', s_ms,'Color',color);hold on; % (left side is stable)
                        CustomMark(x_point,y_point,x,A,color,effective_index);
                    else
                        plot(x_point + delta(1),y_point + delta(2),'o','LineWidth', 1,'MarkerSize', us_ms,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',color);
                        CustomMark(x_point,y_point,x,A,color,effective_index);
                    end
                end

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
                [x_point,y_point] = UVW_to_XY(x);
                % always plot the internal fixed point (eff_index=0)
                CustomMark(x_point,y_point,x,A,color,0); 
            end
        end
    end
    
    marker_legend(effective_index,color);
    
    add_labels(labels);
end

function marker_legend(effective_index,color)
    if (effective_index==0)
        
        % order:saddle, semisource, source (these 6 could also be reversed).

        dark = [0,0,0]+0.5;

        %sink
        y0 = 0.9;
        fixed_point_marker(0,-0.1,y0,color,color,color);
        text([-0.05],[y0],'Sink', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );

        % one zero, one negative
        y0 = 0.85;
        fixed_point_marker(0,-0.1,y0,color,color,dark);
        text([-0.05],[y0],'Semi-sink', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );

        % two zero eigenvalues
        y0 = 0.8;
        fixed_point_marker(0,-0.1,y0,color,dark,dark);
        text([-0.05],[y0],'Neutral', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );

        % saddle
        y0 = 0.75;
        fixed_point_marker(0,-0.1,y0,color,color,[1,1,1]);
        text([-0.05],[y0],'Saddle', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );

        
        
        % one zero, one positive
        y0 = 0.7;
        fixed_point_marker(0,-0.1,y0,color,[1,1,1],dark);
        text([-0.05],[y0],'Semi-source', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );
        
        %source
        y0 = 0.65;
        fixed_point_marker(0,-0.1,y0,color,[1,1,1],[1,1,1]);
        text([-0.05],[y0],'Source', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'fontsize', 14 );

        
    end

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
        type = -1; % zero in the other file
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

% Adapted from:
% - Salman Mashayekh (2020). Custom Marker Plot (https://www.mathworks.com/matlabcentral/fileexchange/39487-custom-marker-plot), MATLAB Central File Exchange. Retrieved July 28, 2020.

function [] = CustomMark(xData,yData,x,A,color,index)
    type = DetermineFixedPointType(x,A);    
    
    color1 = color;
    color2 = [0,0,0]+0.3;
    
    if (type == 4)
        % source: filled circle of white
        color1 = [1,1,1];
        color2 = [1,1,1];
        
    elseif (type == 1)
        % attractor (sink): filled circle of color
        color1 = color;
        color2 = color;
        
    elseif (type == 3)
        % saddle: color + white
        color1 = color;
        color2 = [1,1,1];
        
    elseif (type == -1)
        % two zero eigenvalues
        color1 = [0,0,0]+0.3;
        color2 = [0,0,0]+0.3;
        
    elseif (type == 0)
        % one zero eigenvalue, with stable eigenvalue
        color1 = color;
        color2 = [0,0,0]+0.3;
        
    elseif (type == 2)
        % one zero eigenvalue, with unstable eigenvalue
        color1 = color;
        color2 = [1,1,1];
    end

    if (index==0)

        rotate = 0;
        
        if (x(1) == 1)
            rotate = pi/2;
        elseif (x(2) == 1)
            rotate = -pi/3 - pi/2;
        elseif (x(3) == 1)
            rotate = pi/3 + pi/2;
        else
            if (x(1) == 0)
                rotate = pi/2;
            elseif (x(2) == 0)
                rotate = pi/6 + pi;
            elseif (x(3) == 0)
                rotate = -pi/6 - pi;
            end
        end
        
        fixed_point_marker(rotate,xData,yData,color,color1,color2);
    end
end

function [] = fixed_point_marker(rotate,xData,yData,color,color1,color2)
    r = 1/55;
    markerSize = 1;

    lw = 1.5;
    markerEdgeColor = color;
    markerFaceColor = color;

    % color-filled circle, then gray circle
    for i = 1:2
        if i == 2
            r = 1/70;%slightly smaller
            markerFaceColor = color2;
            phi2 = rotate:0.01:(pi+rotate);
            markerDataX = r*cos(phi2);
            markerDataY = r*sin(phi2);
        else
            markerFaceColor = color1;
            phi1 = 0:0.01:(2*pi);
            markerDataX = r*cos(phi1);
            markerDataY = r*sin(phi1);
        end

        xData = reshape(xData,length(xData),1) ;
        yData = reshape(yData,length(yData),1) ;
        markerDataX = markerSize * reshape(markerDataX,1,length(markerDataX)) ;
        markerDataY = markerSize * reshape(markerDataY,1,length(markerDataY)) ;

        vertX = repmat(markerDataX,length(xData),1) ; vertX = vertX(:) ;
        vertY = repmat(markerDataY,length(yData),1) ; vertY = vertY(:) ;

        vertX = repmat(xData,length(markerDataX),1) + vertX ;
        vertY = repmat(yData,length(markerDataY),1) + vertY ;
        faces = 0:length(xData):length(xData)*(length(markerDataY)-1) ;
        faces = repmat(faces,length(xData),1) ;
        faces = repmat((1:length(xData))',1,length(markerDataY)) + faces ;

        patchHndl = patch('Faces',faces,'Vertices',[vertX vertY]);

        if (i>1)
            set(patchHndl,'FaceColor',markerFaceColor,'EdgeColor','none') ;
        else
            set(patchHndl,'FaceColor',markerFaceColor,'EdgeColor',markerEdgeColor,'LineWidth',lw);
        end
    end
end



