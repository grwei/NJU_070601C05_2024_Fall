%% ModWaveNumAnalysis.m
% Description: create figure of the results of modified wave number analysis
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 16, 2024
% Last modified:
%

%% class def

classdef ModWaveNumAnalysis < handle
    properties (Access=private)
        scheme_disp_name_
        mod_wave_func_handle_
    end

    methods (Access=public)
        function obj = ModWaveNumAnalysis(mod_wave_func_handle, scheme_disp_name_)
            arguments
                mod_wave_func_handle = {
                    @(x) 4 * sin(sqrt(x) / 2).^2;
                    @(x) 6 - 18 ./ (3 + 2 * tan(sqrt(x) / 2).*2);
                }
                scheme_disp_name_ = {
                    "$D_2: \; y = 4 \sin^2 (\sqrt{x} / 2)$";
                    "$D_4: \; y = 6 - 18 / (3 + 2 \tan^2 (\sqrt{x} / 2))$";
                }
            end

            obj.mod_wave_func_handle_ = mod_wave_func_handle;
            obj.scheme_disp_name_ = scheme_disp_name_;
        end
        
        function obj = set_scheme_disp_name(obj, scheme_disp_name)
            obj.scheme_disp_name_ = scheme_disp_name;
        end
        
        function obj = set_mod_wave_func_handle(obj, mod_wave_func_handle)
            obj.mod_wave_func_handle_ = mod_wave_func_handle;
        end

        function scheme_disp_name = get_scheme_disp_name(obj)
            scheme_disp_name = obj.scheme_disp_name_; 
        end

        function mod_wave_func_handle = get_mod_wave_func_handle(obj)
            mod_wave_func_handle = obj.mod_wave_func_handle_;
        end
        
        function obj = create_fig(obj, flag_save)
            arguments
                obj 
                flag_save = false;
            end
            t_fig = figure(Name="fig_5_1_mod_wave_num_results");

            % set figure size
            UNIT_ORIGINAL = t_fig.Units;
            t_fig.Units = "centimeters";
            t_fig.Position = [3, 3, 12, 12];
            t_fig.Units = UNIT_ORIGINAL;

            % create figure
            t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="compact");
            t_axes = nexttile(t_TCL, 1);
            
            x = linspace(0, pi^2, 10000);
            hold(t_axes,"on" )
            for i = 1:length(obj.mod_wave_func_handle_)
                plot(x, obj.mod_wave_func_handle_{i}(x), LineWidth=1.5, DisplayName=obj.scheme_disp_name_{i})
            end
            plot([0, pi^2], [0, pi^2], "--k", LineWidth=1.5, DisplayName="exact: $y = x$")
            
            set(t_axes, FontName="Times New Roman", FontSize=10.5, XTick=(pi/4:pi/4:pi).^2, XTickLabel=["$(\pi/4)^2$", "$(\pi/2)^2$", "$(3\pi/4)^2$", "$\pi^2$"], Box="on", TickLabelInterpreter="latex", XLimitMethod="tight")
            set(t_axes, YTick=t_axes.XTick, YTickLabel=t_axes.XTickLabel, YTickLabelRotation=90)
            grid(t_axes, "on")
            legend(t_axes, Location="northwest", Interpreter="latex", FontSize=10.5, Box="off")
            xlabel(t_axes, "$(k \Delta)^2$", Interpreter="latex", FontSize=10.5);
            ylabel(t_axes, "$(k' \Delta)^2$", Interpreter="latex", FontSize=10.5);
            
            if flag_save
                print(t_fig, t_fig.Name + ".svg", "-dsvg")
            end
        end
    end
end
