function [] = add_gridlines(gridlines,varargin)
    
    p = inputParser;
    color = [0,0,0];
    
    % read in optional parameters    
    [nParams] = length(varargin);
    for param = 1:1:(nParams/2)
        ind = (param-1)*2 + 1;        
        if strcmp(varargin{ind}, 'Color')
            color=varargin{ind+1};
        end
    end
    

    step = 1/gridlines;
    for lattice_delta = step:step:(1-step)
        for i = [1,2]
            for d = [1,2]
                j = mod(i-1 + d,3)+1;
                if (i < j)
                    
                    line_x = zeros(1,3);

                    %% line data for offset line for two-by-two sub-game:                    
                    line_x(1,:) = zeros(1,3);
                    line_x(1,i) = lattice_delta;
                    line_x(1,j) = 1-lattice_delta;

                    if (i == 1) && (j == 3)
                        line_x(2,:) = zeros(1,3) + 1-lattice_delta;
                        line_x(2,i) = lattice_delta;
                        line_x(2,j) = 0;
                        
                        
                        
                    else
                        line_x(2,:) = zeros(1,3) + lattice_delta;
                        line_x(2,i) = 0;
                        line_x(2,j) = 1-lattice_delta;
                        
                
                    end

                    [x_line,y_line] = UVW_to_XY(line_x);
                    plot(x_line,y_line,'-', 'LineWidth', 1,'Color',color);hold on; 
                end
            end
        end
    end

    add_labels({'','',''});

end