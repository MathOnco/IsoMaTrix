function [] = HAL_line_plot(filepath)

    %% unspecified parameters    
%     if (~exist('color','var') || isempty(color))
%         color = [0,0,0]; % black if none specified
%     end
%     
    h = gcf;
    figure_number=h.Number;
    figure(figure_number); hold on;
    
    if (~exist('filepath','var') || isempty(filepath))
        filepath  = 'HALMatrixGame/HALMatrix-output/HAL_trajectory.csv';
    end
    
    %% read in data
    data = dlmread(filepath,',',1,0);
    x_over_time = data(:,1:3);
    x_over_time = x_over_time./sum(x_over_time(1,:));
    [tF,~] = size(x_over_time);
    t = 0:1:tF-1;
    
    plot(t,x_over_time(:,1),'-','LineWidth',2); hold on;
    plot(t,x_over_time(:,2),'-','LineWidth',2);
    plot(t,x_over_time(:,3),'-','LineWidth',2);
    
    
    %% get hal payoff, make this into payoff
    A = HAL_payoff();
    
    x0 = x_over_time(1,:)
    
    [t, x]=ode45(@(t,n)replicator(t,n,A), [0 tF-1], x0);
    plot(t,x(:,1),':','LineWidth',2); hold on;
    plot(t,x(:,2),':','LineWidth',2);
    plot(t,x(:,3),':','LineWidth',2);
    
    
%     h=legend(char(labels(1)),char(labels(2)),char(labels(3)));
%     set(h,'Location','Best');
        
    h=xlabel('time, t'); set(h,'Interpreter','Latex', 'FontSize', 25);
    h=ylabel('$x_i(t)$'); set(h,'Interpreter','Latex', 'FontSize', 25);
    set(gcf,'color','w');
    
    
    

end