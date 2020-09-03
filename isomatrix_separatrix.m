function[area_vec]=isomatrix_separatrix(A,varargin)

    p = inputParser;
    
    %% set up default values for optional parameters: ('Color' and 'Labels')
    color = [0,0,0]; % this is the line color
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
    
    % define the colors:
    blue = [0.2188,0.4531,0.6914];
    red = [0.7734,0.2188,0.1719];
    green = [0.3086,0.6211,0.2227];
    
    % begin with 100% red, then subtract as we add more basins:
    area_vec = [1,0,0]; 
    
    % for reference later
    x2_0 = [0,1,0]';
    x1_0 = [1,0,0]';
    x3_0 = [0,0,1]';
    full_area = [x1_0,x2_0,x3_0];

    S={[],[],[],[],[],[],[]};% store the seperatrices
    fp = {[],[],[],[],[],[],[]};% store the associated fixed points
    [xedge,yedge]=checkedge(A);%information of fixedpoint on the edge
    if all(yedge(:)==0)
        % but not without color full:
        % plot the full area:            
        [x_points,y_points] = UVW_to_XY(full_area');
        FillArea(x_points, y_points, red);
        add_labels(labels);
        return
    end

    [xmid,ymid]=checkmid(A);%information of fixedpoint in the middle

    for i=1:3
        if yedge(i)==1
            [D,lambda,V]=hessian(xedge(:,i)',A);% hessian matrix, eigen values and eigen vectors of fixed points on the edge.
            if lambda(1,1)*lambda(2,2)<0   % check if it's saddle point
                if abs(V(i,1))>abs(V(i,2)) % getting the eigen vector which is not along the edge
                    k=1;
                else
                    k=2;
                end                
                
                lambdak=lambda(k,k); % corresponding the eigen value
                Vk=V(:,k);
                
                if Vk(i)<0           % make sure the eigen vector goes towards inside of the triangle
                    Vk=-1*Vk;
                end                

                e=0.01;             % small shift along the seperetrix
                if ymid==1
                   e=min([e,min(abs(xmid-xedge(:,i)))]);
                end
                
                x20=xedge(:,i)+e*Vk;  % new intial condition to generate seperatrix. shifted from the saddle point
                while not(all(x20(:)>0)) % make sure intial condition is valid
                    e=e/2;
                    x20=xedge(:,i)+e*V;
                end
                
                x_overtime=traj(x20',A,lambdak); % generate the seperetrix. eigen value determines you're going forward or backtracking.
                S{i}=[xedge(:,i) x_overtime]; % store the seperatrix staring from the saddle point, ending at xend.
            end
        end

    end


    

    if ymid==1
        [~,lambda,V]=hessian(xmid',A);                   % checking whether it's a saddle point in the middle.
        if imag(lambda(1,1))==0
            if lambda(1,1)*lambda(2,2)<0

                e=0.01;
                e=min([e,min(abs(ymid))/2]);

                x20=xmid+e*V(:,1);
                x_overtime=traj(x20',A,lambda(1,1));
                S{4}=[xmid x_overtime];

                x20=xmid-e*V(:,1);
                x_overtime=traj(x20',A,lambda(1,1));
                S{5}=[xmid x_overtime];

                x20=xmid+e*V(:,2);
                x_overtime=traj(x20',A,lambda(2,2));
                S{6}=[xmid x_overtime];

                x20=xmid-e*V(:,2);
                x_overtime=traj(x20',A,lambda(2,2));
                S{7}=[xmid x_overtime];

            end
        end
    end


    

    for i=1:7
        if ~isempty(S{i})
                x_overtime=S{i};
                x_over=x_overtime(:,end);

                if ymid==1                   % because the seperatrix will end up reaching another fixed point. find where it is using xend.
                    xend=xmid;
                    err=norm(x_over-xend);
                else

                    err=1;
                end

                for j=1:3
                    xvtx=zeros(3,1);      % see if the end point is vetex
                    xvtx(j,1)=1;
                    err2=norm(x_over-xvtx);
                    if err2<err
                       xend=xvtx;
                       err=err2;

                    end

                    err2=norm(x_over-xedge(:,j)); % see if the end point if another fixed point on the edge.
                    if err2<err
                       xend=xedge(:,j);
                       err=err2;

                    end

                end

                S{i}=[x_overtime xend];
                fp{i} = xend; 
        end
    end

    % how many separatrices are there?
    % cases 4 - 7 are for saddle internal:
    % check if only ONE case for 1:3, if two there is an internal non-saddle:
    cases = 0;
    for i = 1:7
        if ~isempty(S{i})
            if (i <= 3)
                cases = cases+1;
            end
        end
    end
    
    % corners of the simplex for easy referencce:
    x20 = [0,1,0]';
    x10 = [1,0,0]';
    x30 = [0,0,1]';

    % plot the full area:            
    [x_points,y_points] = UVW_to_XY(full_area');
    AFULL=FillArea(x_points, y_points, red);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASES = 0 
    if (cases == 0) % this should be one global equil which is on an edge:
        % plot the full area, only 1 basin      


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASES = 1    
    elseif (cases == 1)
        % this is one separatrix between two f.p., each on an edge:
               
        i_focal = 0;
        % find the empty one:
        for i = 1:3
            if ~isempty(S{i})
                i_focal = i;
            end
        end
        
        i = i_focal;
        full_trajectory = S{i}';
        [first,second] = minimizeAtIndex(full_trajectory(1,:), full_trajectory(end,:), 1);

        [type1] = DetermineFixedPointType(first,A);
        [type2] = DetermineFixedPointType(second,A);


        if (max(  min(first), min(second) ) > 0)
            % internal f.p. to edge f.p. w/ single sep.
            
        else
            % ^^ one of these must be a source (4)
            if ((type1 == 4) || (type2 == 4))
                % if first is on the bottom
                if (first(1) == 0)
                    if (second(3)==0)
                        % goes from 1-2 edge to 2-3 edge:
                        bottom_area = [x20,S{i},x20];                            
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;                        
                    else
                        % goes from 1-3 edge to 2-3 edge:
                        bottom_area = [x30,S{i},x30];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    end

                else % neither fixed point is on the bottom: 
                    if (full_trajectory(1,3)==0)
                        % goes from 1-2 edge to 1-3 edge:
                        % left-to-right:
                        % corner -> S -> other corner -> original corner
                        bottom_area = [x20,S{i},x30,x20];  
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    else
                        % right-to-left:
                        % corner -> S -> other corner -> original corner
                        bottom_area = [x30,S{i},x20,x30];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    end    
                end
            end
        end
            

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASES = 2
    elseif (cases == 2)
        % scenario 1:
        % -> internal source, 2 unstable (saddle) edges
        
        % scenario 2:
        % -> there are no internal fixed points, but are 3 unstable edges
        
        % scenario 3: 
        % -> there are no internal fixed points, but are 3 stable edges
        
        % scenario 4:
        % -> internal attractor, 2 stable edges
        
        % scenario 5:
        % -> edge saddle to internal attractor, treat it like 1 separatrix
        
        case_indices = [0,0];

        % which two cases:
        j = 1;
        for i = 1:3
            if ~isempty(S{i})
                case_indices(j) = i;
                j = 2;
            end
        end

        % overlay area below combined two separatraces:
        full_trajectory1 = S{case_indices(1)};
        full_trajectory2 = S{case_indices(2)};
        
        all_fixed_points = [full_trajectory1(:,1)';
                            full_trajectory1(:,end)';
                            full_trajectory2(:,1)';
                            full_trajectory2(:,end)'];
        [bool,bool_i]=AreAnyInternal(all_fixed_points);
        [types] = DetermineFixedPointType(all_fixed_points,A);
        
        if ((~bool) && min(types) < 3)
            % this is scenario 3: 3 stable edges, two separatrices
            % 1 basin of attraction (keep it red)
            
            
        else % could be scenario 1 or 2 or 4 or 5:
            
            if (bool)
                % scenario 1 or 4:
                internal_type = DetermineFixedPointType(all_fixed_points(bool_i,:),A);
                if (internal_type == 1) % attractor
                    % this is scenario 4 (only 1 basin)
                    
                    % or, check if this is scenario 5:
                    % dealing only with full_trajectory1;
                    
                    min1 = max(min(full_trajectory1(:,1)), min(full_trajectory1(:,end)));
                    min2 = max(min(full_trajectory2(:,1)), min(full_trajectory2(:,end)));
                    
                    full_trajectory = full_trajectory1;
                    
                    if min1 > min2
                        % trajectory 1 has the internal pt:
                        full_trajectory = full_trajectory2;
                    end
                    
                    % ensure 1st starts on edge:
                    if min(full_trajectory(:,1))>0
                        full_trajectory = fliplr(full_trajectory);
                    end
                    
                    [first,second] = minimizeAtIndex(full_trajectory(:,1)', full_trajectory(:,end)', 1);
                                        
                    % if first is on the bottom
                    if (first(1) == 0)
                        if (second(3)==0)
                            % goes from 1-2 edge to 2-3 edge:
                            bottom_area = [x20,full_trajectory,x20];                            
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;                        
                        else
                            % goes from 1-3 edge to 2-3 edge:
                            bottom_area = [x30,full_trajectory,x30];
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;
                        end

                    else % neither fixed point is on the bottom: 
                        if (full_trajectory(1,3)==0)
                            % goes from 1-2 edge to 1-3 edge:
                            % left-to-right:
                            % corner -> S -> other corner -> original corner
                            bottom_area = [x20,full_trajectory,x30,x20];  
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;
                        else
                            % right-to-left:
                            % corner -> S -> other corner -> original corner
                            bottom_area = [x30,full_trajectory,x20,x30];
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;
                        end    
                    end
                    
                    
                    
                    
                    
                else
                    % scenario 1:
                    % -> there is an internal fixed point which is not a saddle:
                    % --> 2 basins of attraction


                    % ensure 1st starts on edge:
                    if min(full_trajectory1(:,1))>0
                        full_trajectory1 = fliplr(full_trajectory1);
                    end

                    % ensure 2nd starts in center:
                    if min(full_trajectory2(:,1))==0
                        full_trajectory2 = fliplr(full_trajectory2);
                    end


                    [first,second] = minimizeAtIndex(full_trajectory1(:,1)', full_trajectory2(:,end)', 1);

                    if (first(1) == 0)
                        % one f.p. on the 2-3 edge:
                        if (second(3)==0)
                            % one f.p. on the 1-2 edge:
                            bottom_area = [x20,full_trajectory1,full_trajectory2,x20];
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;
                        else
                            % one f.p. on the 1-3 edge:
                            bottom_area = [x30,full_trajectory1,full_trajectory2,x30];
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue);
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;
                        end
                    else
                        % no f.p. on the bottom:
                        % this trajectory goes from 1-2 to 1-3 edge:
                        bottom_area = [x10,full_trajectory1,full_trajectory2,x10];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    end

                end


            else
                % scenario 2:
                % -> there are no internal fixed points, but are 3 unstable edges
                % --> 3 basins of attraction

                % FIRST SEPARATRIX:
                [first,second] = minimizeAtIndex(full_trajectory1(:,1)', full_trajectory1(:,end)', 1);

                if (first(1) == 0)
                    % one f.p. on the 2-3 edge:
                    if (second(3)==0)
                        % one f.p. on the 1-2 edge:
                        bottom_area = [x20,full_trajectory1,x20];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    else
                        % one f.p. on the 1-3 edge:
                        bottom_area = [x30,full_trajectory1,x30];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        b1=FillArea(x_points, y_points, blue);
                        area_vec(1)=area_vec(1)-b1/AFULL;
                        area_vec(3)=area_vec(3)+b1/AFULL;
                    end
                else
                    % no f.p. on the bottom:
                    % this trajectory goes from 1-2 to 1-3 edge:
                    bottom_area = [x10,full_trajectory1,x10];
                    [x_points,y_points] = UVW_to_XY(bottom_area');
                    b1=FillArea(x_points, y_points, blue);    
                    area_vec(1)=area_vec(1)-b1/AFULL;
                    area_vec(3)=area_vec(3)+b1/AFULL;
                end

                % SECOND SEPARATRIX:
                [first,second] = minimizeAtIndex(full_trajectory2(:,1)', full_trajectory2(:,end)', 1);

                if (first(1) == 0)
                    % one f.p. on the 2-3 edge:
                    if (second(3)==0)
                        % one f.p. on the 1-2 edge:
                        bottom_area = [x20,full_trajectory2,x20];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        g1=FillArea(x_points, y_points, green);
                        area_vec(1)=area_vec(1)-g1/AFULL;
                        area_vec(2)=area_vec(2)+g1/AFULL;
                    else
                        % one f.p. on the 1-3 edge:
                        bottom_area = [x30,full_trajectory2,x30];
                        [x_points,y_points] = UVW_to_XY(bottom_area');
                        g1=FillArea(x_points, y_points, green);
                        area_vec(1)=area_vec(1)-g1/AFULL;
                        area_vec(2)=area_vec(2)+g1/AFULL;
                    end
                else
                    % no f.p. on the bottom:
                    % this trajectory goes from 1-2 to 1-3 edge:
                    bottom_area = [x10,full_trajectory2,x10];
                    [x_points,y_points] = UVW_to_XY(bottom_area');
                    g1=FillArea(x_points, y_points, green); 
                    area_vec(1)=area_vec(1)-g1/AFULL;
                    area_vec(2)=area_vec(2)+g1/AFULL;
                end
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASES = 3
    elseif (cases == 3)

        % scenario 1:
        % -> stable internal, with stable f.p.'s on every edge
        % ---> this case is all one basin (red)

        % scenario 2:
        % -> 3 edge fixed points and 1 center fixed point which is a source

        [type] = DetermineFixedPointType(xmid',A);
        if (type == 4)

            % make them all start on edge:
            if min(S{1}(:,1) > 0)
                S{1} = fliplr(S{1});
            end

            if min(S{2}(:,1) > 0)
                S{2} = fliplr(S{2});
            end

            if min(S{3}(:,1) > 0)
                S{3} = fliplr(S{3});
            end

            case_indices = [0,0,0];

            for i = 1:3

                % check first entry:
                x0 = S{i}(:,1);

                if x0(3) == 0
                    % this one begins on 1-2:
                    case_indices(1) = i;
                end

                if x0(2) == 0
                    % this one begins on 1-3:
                    case_indices(2) = i;
                end

                if x0(1) == 0
                    % this one begins on 2-3:
                    case_indices(3) = i;
                end

            end

            % color blue corner:
            % ---> % use traj 1, with flipped traj 3, and x20
            second_area = [x20,S{case_indices(1)},fliplr(S{case_indices(3)}),x20];
            [x_points,y_points] = UVW_to_XY(second_area');
            b1=FillArea(x_points, y_points, blue);
            area_vec(1)=area_vec(1)-b1/AFULL;
            area_vec(3)=area_vec(3)+b1/AFULL;

            % this basin goes to x20 ^^


            % color green corner:
            % ---> % use traj 2, with flipped traj 3, and x20
            third_area = [x30,S{case_indices(2)},fliplr(S{case_indices(3)}),x30];
            [x_points,y_points] = UVW_to_XY(third_area');
            g1=FillArea(x_points, y_points, green);
            area_vec(1)=area_vec(1)-g1/AFULL;
            area_vec(2)=area_vec(2)+g1/AFULL;

            % this basin goes to x30 ^^

        end
    end

    % final four cases stem out from internal f.p.
    % --> 4 is paired w/ 5; 6 is paired with 7:
    for i = [4, 6]
        if ~isempty(fp{i})
            
            % find trajectories:
            
            % make sure they start on the end:
            % make them all start on edge:
            if min(S{i}(:,1) > 0)
                S{i} = fliplr(S{i});
            end

            if min(S{i+1}(:,1) > 0)
                S{i+1} = fliplr(S{i+1});
            end
            
            all_fixed_points = [S{i}(:,1)';
                                S{i}(:,end)';
                                S{i+1}(:,1)';
                                S{i+1}(:,end)'];
                                
            types = DetermineFixedPointType(all_fixed_points,A);
            
            if ((types(2) == 3) && (types(2) == 3))
                % internal saddle:
                
                % color it blue if it goes UNSTABLE -> saddle -> UNSTABLE
                % i think 4 or 0?
                
                if ((types(1) == 0) || (types(1) == 4))
                    if ((types(3) == 0) || (types(3) == 4))
                        
                        % color blue below both seps:
                        [first,second] = minimizeAtIndex(S{i}(:,1)', S{i+1}(:,1)', 1);
                        
                        % if the midpoint is on left edge, use bottom left corner (x0)
                        if first(1) == 0
                            % one is on 2-3 edge:
                            if (second(3)==0)
                                % the other f.p. is on the 1-2 edge:
                                bottom_area = [x20,S{i},fliplr(S{i+1}),x20];
                                [x_points,y_points] = UVW_to_XY(bottom_area');
                                b1=FillArea(x_points, y_points, blue);
                                area_vec(1)=area_vec(1)-b1/AFULL;
                                area_vec(3)=area_vec(3)+b1/AFULL;                                
                            else
                                % the other f.p. is on the 1-3 edge:
                                bottom_area = [x30,S{i},fliplr(S{i+1}),x30];
                                [x_points,y_points] = UVW_to_XY(bottom_area');
                                b1=FillArea(x_points, y_points, blue);
                                area_vec(1)=area_vec(1)-b1/AFULL;
                                area_vec(3)=area_vec(3)+b1/AFULL;
                            end
                            
                        else
                            % neither is on 2-3 edge:                            
                            % therefore, this trajectory goes from 1-2 to 1-3 edge:
                            bottom_area = [x10,S{i},fliplr(S{i+1}),x10];
                            [x_points,y_points] = UVW_to_XY(bottom_area');
                            b1=FillArea(x_points, y_points, blue); 
                            area_vec(1)=area_vec(1)-b1/AFULL;
                            area_vec(3)=area_vec(3)+b1/AFULL;                            
                        end
                    end
                end
                
            else
                % do nothing (not an internal saddle)
                
            end
        end
    end

    % lastly, draw on separatrices 
    for i = 1:7
        if ~isempty(S{i})
            [x_points,y_points] = UVW_to_XY(S{i}');
            plot(x_points,y_points,linestyle, 'LineWidth', linewidth,'Color',color);hold on;
        end
    end
    
    % TK
%     labels
    add_labels(labels);
    
    
end


function [type] = DetermineFixedPointType(xvec,A)

    [n,~] = size(xvec);
    type = zeros(n,1);

    for i = 1:n
        x=xvec(i,:);

        [~,lambda,~]=hessian(x,A);

        lambda1 = real(lambda(1,1));
        lambda2 = real(lambda(2,2));

        if ((lambda1== 0) && (lambda2==0))
            % indeterminant
            type(i) = 0;
        elseif (((lambda1== 0) && (lambda2>0)) || ((lambda1>0) && (lambda2==0)) )
            % one zero eigenvalue, one negative
            type(i) = 0;
        elseif (((lambda1== 0) && (lambda2<0)) || ((lambda1< 0) && (lambda2==0)) )
            % one zero eigenvalue, one positive (UNSTABLE)
            type(i) = 2;
        elseif lambda1*lambda2<0   
            % this is a saddle point
            type(i) = 3;
        elseif ((lambda1> 0) && (lambda2>0))
            % this is a source
            type(i) = 4;
        else
            % this is an attractor
            type(i) = 1;
        end
    end
end

function [first,second] = minimizeAtIndex(temp_first, temp_second, i)
    first = temp_first;
    second= temp_second;
    % switch their position if one has zero:
    if (temp_first(i) > temp_second(i))
        first = temp_second;
        second= temp_first;
    end
end

function[x,y]=checkedge(B)
    x=zeros(3,3); % each column stores the fixed point on the edge.
    y=zeros(1,3); % whether the fixed point is inside the triangle
    for i=1:3
        x1=1;
        if x1==i
            x1=x1+1;
        end
        x2=x1+1;
        if x2==i
            x2=x2+1;
        end

        if (B(x1,x1)-B(x2,x1))*(B(x1,x2)-B(x2,x2))<0
            y(i)=1;
            x(x1,i)=(B(x2,x2)-B(x1,x2))/(B(x2,x2)-B(x1,x2)+B(x1,x1)-B(x2,x1));
            x(x2,i)=(B(x1,x1)-B(x2,x1))/(B(x2,x2)-B(x1,x2)+B(x1,x1)-B(x2,x1));
        end
    end
end

% check internal f.p.
function [x,y]=checkmid(B)
    x=0;  %fixed point in the middle
    y=0;  % whether the fixed point is inside the triangle

    A=[B(1,:)-B(2,:);
        B(1,:)-B(3,:);
        1 1 1];

    if rank(A)==3
        mid=A\[0;0;1];
        if all(mid(:)>0)
            y=1;
            x=mid;
        end
    end
end

% trajectory along separatrix
function[x_overtime]=traj(x0,B,lambda)
    q=eye(3); % mutation matrix
    w=[1 1 1];% selection pressure, just one
    dt=0.01;  % step side
    step=100/dt;  % number of steps
    if lambda<0   % check if forward or backward
        dt=-1*dt;
    end
    x_overtime=zeros(3,step);
    dif=abs([B(1,:)-B(2,:) B(1,:)-B(3,:)]); % see if we need to amplify the matrix or not
    m=min(dif(dif>0));  % the minimum positive difference between row elements
    amp=1;              % amplification
    if m>0
        amp=0.5/m;
    end

    B=B*amp;            

    for i=1:step
        x0=repli2(x0,B,q,w,dt); % rk45 method of 1 step.
        x_overtime(:,i)=x0; 
    end
end

function [a] = FillArea(x_points, y_points, color)
h = fill(x_points,y_points,color);
h.FaceColor = color;
h.LineStyle = 'none';
a = polyarea(x_points,y_points);
end

% determine if any fixed points are internal
function [bool,bool_i] = AreAnyInternal(fixed_points)
    bool_i = -1;
    [n,~] = size(fixed_points);
    bool = false;
    for i = 1:n
        if min(fixed_points(i,:)>0)
            bool = true;
            bool_i = i;
        end
    end
end


function x2=repli2(x,payoff,q,w,dt)
    x1=x;
    f1=1-w+w.*(payoff*x1')';
    fave=x1*f1';
    k1=f1.*x1*q-fave*x1;
    
    
    x2=x1+dt/2*k1;
    f2=1-w+w.*(payoff*x2')';
    fave=x2*f2';
    k2=f2.*x2*q-fave*x2;

    x3=x1+dt/2*k2;
    f3=1-w+w.*(payoff*x3')';
    fave=x3*f3';
    k3=f3.*x3*q-fave*x3;
    
    
    x4=x1+dt*k3;
    f4=1-w+w.*(payoff*x4')';
    fave=x4*f4';
    k4=f4.*x4*q-fave*x4;
    
    
    x2=x+dt/6*(k1+2*k2+2*k3+k4);
    if sum(x2)~=1
        x2=x2/sum(x2);
    end
end

