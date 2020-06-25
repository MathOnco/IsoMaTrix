function [] = line_plot(A,x0,tF,varargin)

    %% set up default values for optional parameters: ('Labels')
    labels = {'1','2','3'};
    [nn,~] = size(A);
    assert(3==nn,'Please provide a 3 by 3 matrix.')
    p = inputParser;
    
    % read in optional parameters    
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Labels')
            labels=varargin{ind+1};
        end
    end

    
    

    [t, x]=ode45(@(t,n)replicator(t,n,A), [0 tF], x0);
    plot(t,x(:,1),'LineWidth',2); hold on;
    plot(t,x(:,2),'LineWidth',2);
    plot(t,x(:,3),'LineWidth',2);
    
    
    h=legend(char(labels(1)),char(labels(2)),char(labels(3)));
    set(h,'Location','Best');
    h=xlabel('time, t'); set(h,'Interpreter','Latex', 'FontSize', 25);
    h=ylabel('$x_i(t)$'); set(h,'Interpreter','Latex', 'FontSize', 25);
    set(gcf,'color','w');ylim([0 1]);xlim([0 tF]);
end