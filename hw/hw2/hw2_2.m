%% hw2_2.m
% Description: Solving the 1-D advection problem: c_t + u c_x = 0
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 21, 2024
% Last modified: Oct. 24, 2024
%

clear; clc; close all

%% Problem 2(b)

[f_list, x_grid, t_list] = hw2_2b("time_forward_Euler_space_upwind", true);

%% Problem 2(c)

hw2_2b("time_forward_Euler_space_central_diff", true); % Problem 2(c)

%% Function: Problem 2(b, c)

function [f_list, x_grid, t_list] = hw2_2b(scheme_name, flag_save)
    arguments (Input)
        scheme_name = "time_forward_Euler_space_upwind";
        flag_save = false;
    end
    
    solver = AdvecProb2D();

    %%% Solver params (1-D Advection problem)

    velocity = [1, ];
    t_start = 0;
    x_range = {[0, 4], };
    init_func = @(x) double(x > 0 & ~(x > .5));
    bndry_func = @(x, t) double(x > x_range{1}(1) & x < x_range{1}(2));
    delta_x = 0.05;
    t_query = [0, 1, 2, 3];

    %%% Solving

    Cr_list = [.5, .99, 1, 1.01]; % Cr = u delta_t / delta_x
    f_list = cell(length(Cr_list), 1);
    x_grid = f_list;
    t_list = x_grid;
    for test_idx = 1:length(Cr_list)
        Cr = Cr_list(test_idx);
        delta_t =  Cr * delta_x / velocity(1);
        params_struct = AdvecProb2D.prepare_params(velocity, t_start, x_range, init_func, bndry_func, delta_t, delta_x, t_query);
        solver.reset(params_struct);
        solver.solve(scheme_name);
        [f_list{test_idx}, x_grid{test_idx}, t_list{test_idx}, ~] = solver.get_solution();
    end

    %%% Create figure
    switch scheme_name
        case "time_forward_Euler_space_upwind"
            fig_name_prefix = "fig_2_upwind";
            scheme_name_short = "upwind";
        case "time_forward_Euler_space_central_diff"
            fig_name_prefix = "fig_2_centdif";
            scheme_name_short = "central diff.";
        otherwise
            fig_name_prefix = "fig_2_undef";
            scheme_name_short = "unknown scheme";
    end

    t_fig = figure(Name=sprintf("%s_Cr_compare_1", fig_name_prefix));

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 16, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, floor(length(Cr_list) / 2), 2, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$", Interpreter="latex", FontSize=10.5);
    ylabel(t_TCL, "(numerical, " + scheme_name_short + ") $f$ (subjected to $f_{t} + u f_x = 0, \; f|_{t = 0} = f^0, \; f|_{x \in \partial \Omega} = 0$)", Interpreter="latex", FontSize=10.5);
    for test_idx = 1:length(Cr_list)
        t_axes = nexttile(t_TCL, test_idx);
        hold(t_axes, "on");
        for t_query_idx = 1:length(t_list{test_idx})
            plot(t_axes, x_grid{test_idx}{1}, f_list{test_idx}{t_query_idx}, LineWidth=1.0, DisplayName=sprintf("$t = %.2f$", t_list{test_idx}(t_query_idx)))
        end
        set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
            Tag=sprintf("(%c) $\\mathrm{Cr} = %g$", char('a' - 1 + test_idx), Cr_list(test_idx)));
        if scheme_name == "time_forward_Euler_space_upwind"
            t_axes.XLim = x_range{1};
            t_axes.YLim = [0, 2];
        end
        if scheme_name == "time_forward_Euler_space_upwind" && ~mod(test_idx, 2)
            t_axes.YTickLabel = {};
        end
        if scheme_name == "time_forward_Euler_space_upwind" && test_idx < length(Cr_list) - 1
            t_axes.XTickLabel = {};
        end
        if scheme_name == "time_forward_Euler_space_upwind" 
            legend(t_axes, Location="northwest", Interpreter="latex", FontSize=10.5, Box="off");
        elseif scheme_name == "time_forward_Euler_space_central_diff" 
            legend(t_axes, Location="southeast", Interpreter="latex", FontSize=10.5, Box="off");
        else
            legend(t_axes, Location="best", Interpreter="latex", FontSize=10.5, Box="off");
        end
    end

    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        t_txt_box = annotation(t_fig, "textbox", String="\bf " + t_axes.Tag, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*10, 10.5*1.5], 10.5*10, 10.5*1.5];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_save
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")
    end

    %%% Fig. Cr compare fixed t
    t_fixed = 2;
    t_fig = figure(Name=sprintf("%s_Cr_compare_2", fig_name_prefix));

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 12, 12];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="compact");
    t_axes = nexttile(t_TCL, 1);
    
    test_idx = find(Cr_list == 1, 1);
    if ~isempty(test_idx)
        [~, t_query_idx] = min(t_list{test_idx} - t_fixed, [], ComparisonMethod="abs");
        plot(t_axes, x_grid{test_idx}{1}, f_list{test_idx}{t_query_idx}, "-k", LineWidth=1.0, DisplayName=sprintf("$\\mathrm{Cr} = %g, \\; t = %.2f$", Cr_list(test_idx), t_list{test_idx}(t_query_idx)));
    end
    hold(t_axes, "on");
    for test_idx = 1:length(Cr_list)
        if Cr_list(test_idx) == 1
            continue
        end
        [~, t_query_idx] = min(t_list{test_idx} - t_fixed, [], ComparisonMethod="abs");
        plot(t_axes, x_grid{test_idx}{1}, f_list{test_idx}{t_query_idx}, LineWidth=1.0, DisplayName=sprintf("$\\mathrm{Cr} = %g, \\; t = %.2f$", Cr_list(test_idx), t_list{test_idx}(t_query_idx)));
    end

    set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", XLimitMethod="tight");
    legend(t_axes, Location="best", Interpreter="latex", FontSize=10.5, Box="off");

    xlabel(t_axes, "$x$", Interpreter="latex", FontSize=10.5);
    ylabel(t_axes, "(numerical, " + scheme_name_short + ") $f$", Interpreter="latex", FontSize=10.5);

    if flag_save
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")
    end

end
