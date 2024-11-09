%% hw3_1.m
% Description: Solving the 1-D advection problem: c_t + u c_x = 0
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Nov. 7, 2024
% Last modified: Nov. 9, 2024
%

clear; clc; close all

addpath("..\hw2\")
if ~isfolder(".\fig\")
    mkdir(".\fig\")
end

%% Problem 1(a)

Cr_list = [.50, .99, 1.00, 1.01]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=0, alpha=1);

%% Problem 1(d)

Cr_list = [.50, .75, .975, 1]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=.05, alpha=1);

Cr_list = [.50, .75, .975, 1]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=.05, alpha=.53);

Cr_list = [.50, .25, .57, .60]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=1, alpha=1);

Cr_list = [.50, .25, .57, .60]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=1, alpha=.53);

Cr_list = [.50, .25, .57, .60]; % Cr = u delta_t / delta_x
[f_list, x_grid, t_list] = hw3_1_unit("time_leap_frog_RA_space_central_diff", Cr_list, true, nu=.2, alpha=.53);

%% 首步扰动试验

clc;clear;close all;

scheme_name = "time_leap_frog_RA_space_central_diff";
Cr_list = [.5, .99, 1, 1.01];
delta_x = 0.25;
x_range = {[0, delta_x*4*(2^2)], };
% init_func = @(x) double(x > 0 & ~(x > .5));
init_func = @(x) zeros(size(x));
bndry_func = @(x, t) "periodic";
t_query = [0, 1, 2, 4, 8, 16];
ra_opts.nu = 0;
ra_opts.alpha = 1;
flag_save = true;

%%% 周期 4h

disturb_func_first_step = @(x) 1*(sin(2*pi / 4 / delta_x * x) + cos(2*pi / 4 / delta_x * x));
disturb_func_disp_str = "$f_{\rm{dist}} = 1 \cdot (\sin{(2 \pi x / (4h))} + \cos{(2 \pi x / (4h))})$";
fig_name = "fig_disturb_test_4h";

[f_list, x_grid, t_list] = disturb_test_unit(scheme_name, Cr_list, x_range, init_func, bndry_func, disturb_func_first_step, disturb_func_disp_str, delta_x, t_query, flag_save, fig_name, nu=ra_opts.nu, alpha=ra_opts.alpha);

%%% 周期 2h (Nyquist)

disturb_func_first_step = @(x) 1*(sin(2*pi / 2 / delta_x * x) + cos(2*pi / 2 / delta_x * x));
disturb_func_disp_str = "$f_{\rm{dist}} = 1 \cdot (\sin{(2 \pi x / (2h))} + \cos{(2 \pi x / (2h))})$";
fig_name = "fig_disturb_test_2h";

[f_list, x_grid, t_list] = disturb_test_unit(scheme_name, Cr_list, x_range, init_func, bndry_func, disturb_func_first_step, disturb_func_disp_str, delta_x, t_query, flag_save, fig_name, nu=ra_opts.nu, alpha=ra_opts.alpha);

%%% 周期 8h

disturb_func_first_step = @(x) 1*(sin(2*pi / 8 / delta_x * x) + cos(2*pi / 8 / delta_x * x));
disturb_func_disp_str = "$f_{\rm{dist}} = 1 \cdot (\sin{(2 \pi x / (8h))} + \cos{(2 \pi x / (8h))})$";
fig_name = "fig_disturb_test_8h";

[f_list, x_grid, t_list] = disturb_test_unit(scheme_name, Cr_list, x_range, init_func, bndry_func, disturb_func_first_step, disturb_func_disp_str, delta_x, t_query, flag_save, fig_name, nu=ra_opts.nu, alpha=ra_opts.alpha);

%% Function: Unit

function [f_list, x_grid, t_list] = hw3_1_unit(scheme_name, Cr_list, flag_save, ra_opts)
    arguments (Input)
        scheme_name = "time_leap_frog_RA_space_central_diff";
        Cr_list = [.5, .99, 1, 1.01]; 
        flag_save = false;
        ra_opts.nu = 0;
        ra_opts.alpha = 1;
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

    % Cr_list = [.5, .99, 1, 1.01]; % Cr = u delta_t / delta_x
    f_list = cell(length(Cr_list), 1);
    x_grid = f_list;
    t_list = x_grid;
    for test_idx = 1:length(Cr_list)
        Cr = Cr_list(test_idx);
        delta_t =  Cr * delta_x / velocity(1);
        params_struct = AdvecProb2D.prepare_params(velocity, t_start, x_range, init_func, bndry_func, delta_t, delta_x, t_query);
        solver.reset(params_struct);
        solver.solve(scheme_name, nu=ra_opts.nu, alpha=ra_opts.alpha);
        [f_list{test_idx}, x_grid{test_idx}, t_list{test_idx}, ~] = solver.get_solution();
    end

    %%% Create figure
    switch scheme_name
        case "time_forward_Euler_space_upwind"
            fig_name_prefix = "fig_upwind";
            scheme_name_short = "upwind";
        case "time_forward_Euler_space_central_diff"
            fig_name_prefix = "fig_centdif";
            scheme_name_short = "central diff.";
        case "time_leap_frog_RA_space_central_diff"
            fig_name_prefix = sprintf("fig_LF_RA_CD_nu_%.2g_alpha_%.2g", ra_opts.nu, ra_opts.alpha);
            scheme_name_short = "LF\_RA\_CD";
        otherwise
            fig_name_prefix = "fig_undefined";
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
    if scheme_name == "time_leap_frog_RA_space_central_diff"
        title(t_TCL, {"$c^n_i \leftarrow c_i^n + \alpha d_f, \; c^{n+1}_i \leftarrow c_i^{n+1} - (1 - \alpha) d_f, \; d_f = \frac{\nu}{2} \left( c_{i}^{n+1} - 2c_i^{n} + c_{i}^{n-1} \right)$"; sprintf("$\\nu$ = %.3g, $\\alpha$ = %.3g", ra_opts.nu, ra_opts.alpha)}, Interpreter="latex", FontSize=10.5);
    end

    for test_idx = 1:length(Cr_list)
        t_axes = nexttile(t_TCL, test_idx);
        hold(t_axes, "on");
        for t_query_idx = 1:length(t_list{test_idx})
            plot(t_axes, x_grid{test_idx}{1}, f_list{test_idx}{t_query_idx}, LineWidth=1.0, DisplayName=sprintf("$t = %.2f$", t_list{test_idx}(t_query_idx)))
        end
        set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
            Tag=sprintf("(%c) $\\mathrm{Cr} = %g$", char('a' - 1 + test_idx), Cr_list(test_idx)));
        legend(t_axes, Location="northwest", Interpreter="latex", FontSize=10.5, Box="off");
    end

    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        str_{1} = "\bf " + t_axes.Tag;
        t_txt_box = annotation(t_fig, "textbox", String=str_, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*10, 10.5*1.5], 10.5*10, 10.5*1.5];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_save
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-vector", "-dsvg")
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
    if scheme_name == "time_leap_frog_RA_space_central_diff"
        title(t_TCL, {"$c^n_i \leftarrow c_i^n + \alpha d_f, \; c^{n+1}_i \leftarrow c_i^{n+1} - (1 - \alpha) d_f, \; d_f = \frac{\nu}{2} \left( c_{i}^{n+1} - 2c_i^{n} + c_{i}^{n-1} \right)$"; sprintf("$\\nu$ = %.3g, $\\alpha$ = %.3g", ra_opts.nu, ra_opts.alpha)}, Interpreter="latex", FontSize=10.5);
    end

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
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-vector", "-dsvg")
    end

end

%% Func: 首步扰动试验

function [f_list, x_grid, t_list] = disturb_test_unit(scheme_name, Cr_list, x_range, init_func, bndry_func, disturb_func_first_step, disturb_func_disp_str, delta_x, t_query, flag_save, fig_name, ra_opts)
    arguments (Input)
        scheme_name = "time_leap_frog_RA_space_central_diff";
        Cr_list = [.5, .99, 1, 1.01];
        x_range = {[0, 2], };
        init_func = @(x) 5*double(x > 0 & ~(x > .5));
        bndry_func = @(x, t) "periodic";
        disturb_func_first_step = @(x) sin(pi / 2 / obj.delta_x * x) + cos(pi / 2 / obj.delta_x * x);
        disturb_func_disp_str = "$f_{\rm{dist}} = \sin{(\pi x / (2h))} + \cos{(\pi x / (2h))}$";
        delta_x = 0.05;
        t_query = [0, 1, 2, 3];
        flag_save = false;
        fig_name = "fig_disturb_test";
        ra_opts.nu = 0;
        ra_opts.alpha = 1;
    end
    
    solver = AdvecProb2D();

    %%% Solver params (1-D Advection problem)

    velocity = [1, ];
    t_start = 0;
    % x_range = {[0, 4], };
    % init_func = @(x) double(x > 0 & ~(x > .5));
    % bndry_func = @(x, t) double(x > x_range{1}(1) & x < x_range{1}(2));
    % delta_x = 0.05;
    % t_query = [0, 1, 2, 3];

    %%% Solving

    % Cr_list = [.5, .99, 1, 1.01]; % Cr = u delta_t / delta_x
    f_list = cell(length(Cr_list), 1);
    x_grid = f_list;
    t_list = x_grid;
    for test_idx = 1:length(Cr_list)
        Cr = Cr_list(test_idx);
        delta_t =  Cr * delta_x / velocity(1);
        params_struct = AdvecProb2D.prepare_params(velocity, t_start, x_range, init_func, bndry_func, delta_t, delta_x, t_query);
        solver.reset(params_struct);
        solver.solve(scheme_name, disturb_func_first_step, nu=ra_opts.nu, alpha=ra_opts.alpha);
        [f_list{test_idx}, x_grid{test_idx}, t_list{test_idx}, ~] = solver.get_solution();
    end

    %%% Create figure
    switch scheme_name
        case "time_forward_Euler_space_upwind"
            scheme_name_short = "upwind";
        case "time_forward_Euler_space_central_diff"
            scheme_name_short = "central diff.";
        case "time_leap_frog_RA_space_central_diff"
            scheme_name_short = "LF\_RA\_CD";
        otherwise
            scheme_name_short = "unknown scheme";
    end

    t_fig = figure(Name=fig_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 16, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, floor(length(Cr_list) / 2), 2, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$", Interpreter="latex", FontSize=10.5);
    ylabel(t_TCL, "(numerical, " + scheme_name_short + ") $f$ (subjected to $f_{t} + u f_x = 0, \; f|_{t = 0} = f^0, \; f|_{x \in \partial \Omega} = 0$)", Interpreter="latex", FontSize=10.5);
    if scheme_name == "time_leap_frog_RA_space_central_diff"
        title(t_TCL, {"$c^n_i \leftarrow c_i^n + \alpha d_f, \; c^{n+1}_i \leftarrow c_i^{n+1} - (1 - \alpha) d_f, \; d_f = \frac{\nu}{2} \left( c_{i}^{n+1} - 2c_i^{n} + c_{i}^{n-1} \right)$"; ...
                      sprintf("$\\nu$ = %.3g, $\\alpha$ = %.3g, $\\quad$ %s", ra_opts.nu, ra_opts.alpha, disturb_func_disp_str) ...
                      }, Interpreter="latex", FontSize=10.5);
    end

    for test_idx = 1:length(Cr_list)
        t_axes = nexttile(t_TCL, test_idx);
        hold(t_axes, "on");
        for t_query_idx = 1:length(t_list{test_idx})
            plot(t_axes, x_grid{test_idx}{1}, f_list{test_idx}{t_query_idx}, LineWidth=1.0, DisplayName=sprintf("$t = %.2f$", t_list{test_idx}(t_query_idx)))
        end
        set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
            Tag=sprintf("(%c) $\\mathrm{Cr} = %g$", char('a' - 1 + test_idx), Cr_list(test_idx)));
        legend(t_axes, Location="northwest", Interpreter="latex", FontSize=10.5, Box="off");
    end

    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        str_{1} = "\bf " + t_axes.Tag;
        t_txt_box = annotation(t_fig, "textbox", String=str_, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*10, 10.5*1.5], 10.5*10, 10.5*1.5];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_save
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-vector", "-dsvg")
    end
end
