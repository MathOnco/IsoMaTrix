function [] = add_labels(label_string)

    if (~exist('label_string','var') || isempty(label_string))
        label_string = {'1','2','3'};
    end
    
    L = 1;
    H = (L/2)*tan(60*pi/180);
    %% label 1
    delta = 0.1; % was 0.07
    text([L/2,-delta,L+delta],[H+delta,-tan(pi/6)*delta,-tan(pi/6)*delta],label_string, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 36 );
    
    %% plot triangle boundaries:
    lw = 5;
    H = (1/2)*tan(60*pi/180);
    plot([0,1/2],[0,H],'-k','LineWidth',lw);
    plot([1/2,1],[H,0],'-k','LineWidth',lw);
    plot([0,1],[0,0],'-k','LineWidth',lw);
    plot([0,1/2,1],[0,H,0],'.k','MarkerSize',15);

    del = (1 - H)/2; % 0.067
    mul = 1;
    xl = xlim;
    
    xlim([min(-mul*del,xl(1)),max(L + mul*del,xl(2))]);
    
    % set ylim based on xlim, if they've been previously set elsewhere
    xl = xlim;
    del = abs(xl(1));
    mul = max(del/((1 - H)/2),2);
    
    del = (1 - H)/2; % 0.067
    ylim([-mul*del-del,H+mul*del]);
    xlim([min(-mul*del,xl(1)),max(L + mul*del,xl(2))]);

    %% clean up triangle
    set(gca,'visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on'); % leave on title
    box off;
    set(gcf,'color','w');
    pos=get(gcf,'Position');
    width = 600;
    set(gcf,'Position',[pos(1),pos(2),width/sin(pi/3),width]);
end