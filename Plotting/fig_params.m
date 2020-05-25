classdef fig_params  
    properties            
        %Global figure options
        font_size = 16;
        font_name = 'Arial';
        font_weight = 'normal';
        units = 'centimeters';
        line_width = 1;            
        default_color = [0.5 0.5 0.5]; %for any default color plots. 
        
    end

    methods
        function FormatAxes(obj,ax_handle)
            set(ax_handle,'fontsize',obj.font_size,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'linewidth',obj.line_width,...
                'fontsize',obj.font_size,...
                'box','off');
        end %end function
        
        function SetTitle(obj,ax_handle,title_str)
            title(ax_handle,title_str,...
                'fontname',obj.font_name,...
                'fontweight',obj.font_weight,...
                'fontsize',obj.font_size)
        end%end function
        
        function FigureSizing(obj,handles,ax_position,fig_position)
            for i = 1:numel(handles)
                set(0, 'currentfigure', handles(i)); 
                set(gca,'units','centimeters',...
                    'position',ax_position)
                if ~isempty(fig_position)
                    set(handles(i),'units','centimeters',...
                    'position',fig_position)
                end
            end
        end%end function
        
        function SaveFigs(obj,handles,filetype,name,savedir,fancyflag)                
            saveCurFigs(handles,filetype,name,savedir,fancyflag)
            close(handles)
        end%end function
        
        function col = GenDefaultColorArray(obj,array_size)
            col = cell(1,array_size);
            for i = 1:array_size
                col{i} = obj.default_color;
            end
        end %end function
       
       
        
    end

end











