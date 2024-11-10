%% hw4_1.m
% Description: Solving the 1-D advection problem: c_t + u c_x = 0
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Nov. 8, 2024
% Last modified: Nov. 10, 2024
%

clear; clc; close all

addpath("..\hw2\")
if ~isfolder(".\fig\")
    mkdir(".\fig\")
end

flag_save_fig = true;

%% Problem 1(a)

scheme_name = "time_SI_space_CD";
theta_list = [0, .5, .75, 1]; % delta_t = r_num * delta_x^2 / diffu_coeff;
delta_x = .05;
r_num = .45;
init_func = @(x) double(abs(x) < .3);
init_func_disp_str = "$p(x) = 1, \;\; (|x| < 0.3)$";
bndry_func_type = "Neumann";
bndry_func = @(x, t) zeros(size(x));
bndry_func_disp_str = "$q(x) = 0$";
fig_name = "fig_1_Neumann";
[f_list, x_grid, t_list] = hw4_1_unit(scheme_name, theta_list, delta_x, r_num, bndry_func_type, bndry_func, bndry_func_disp_str, init_func, init_func_disp_str, flag_save_fig, fig_name);

%% Problem 1(b)

%%% p(x) = sin(2 \pi x)

scheme_name = "time_SI_space_CD";
theta_list = [0, .5, .75, 1]; % delta_t = r_num * delta_x^2 / diffu_coeff;
delta_x = .05;
r_num = .45;
init_func = @(x) sin(2*pi*x);
init_func_disp_str = "$p(x) = \sin{(2 \pi x)}$";
bndry_func_type = "Dirichlet";
bndry_func = @(x, t) zeros(size(x));
bndry_func_disp_str = "$q(x) = 0$";
% flag_save_fig = true;
fig_name = "fig_2_Dirichlet_1";
[f_list, x_grid, t_list] = hw4_1_unit(scheme_name, theta_list, delta_x, r_num, bndry_func_type, bndry_func, bndry_func_disp_str, init_func, init_func_disp_str, flag_save_fig, fig_name);

%%% p(x) = sin(1 \pi x)

init_func = @(x) sin(1*pi*x);
init_func_disp_str = "$p(x) = \sin{(1 \pi x)}$";
fig_name = "fig_2_Dirichlet_2";
[f_list, x_grid, t_list] = hw4_1_unit(scheme_name, theta_list, delta_x, r_num, bndry_func_type, bndry_func, bndry_func_disp_str, init_func, init_func_disp_str, flag_save_fig, fig_name);

%%% p(x) = sin(4 \pi x)

init_func = @(x) sin(4*pi*x);
init_func_disp_str = "$p(x) = \sin{(4 \pi x)}$";
fig_name = "fig_2_Dirichlet_3";
[f_list, x_grid, t_list] = hw4_1_unit(scheme_name, theta_list, delta_x, r_num, bndry_func_type, bndry_func, bndry_func_disp_str, init_func, init_func_disp_str, flag_save_fig, fig_name);

%% SI_CD 格式的 von Neumann 稳定度分析

clc; clear;

r_vec = 0:.001:2;
theta_vec = 0:.001:1;
[r_grid, theta_grid] = ndgrid(r_vec, theta_vec);
G_max_grid = arrayfun(@(r, theta) G_abs_max(r, theta), r_grid, theta_grid);

%% create figure

t_fig = figure(Name="fig_3_SI_CD_stability");
    
% set figure size
UNIT_ORIGINAL = t_fig.Units;
t_fig.Units = "centimeters";
t_fig.Position = [3, 3, 16, 16];
t_fig.Units = UNIT_ORIGINAL;
t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="tight");
t_axes = nexttile(t_TCL, 1);

[C, h] = contour(t_axes, theta_vec, r_vec, G_max_grid, LevelList=[1.0001, 1.0001], LineWidth=2.0);
clabel(C, h, FontSize=8, Interpreter="latex");
hold on
[C, h] = contour(t_axes, theta_vec, r_vec, G_max_grid, LevelList=0:.25:2, LineWidth=1.0, DisplayName="$$\max_{|\xi h| \le \pi} |G|$$");
hold off
clabel(C, h, [.5, 1.0, 1.5], FontSize=8, Interpreter="latex");

set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", Tag="$|G|_{\mathrm{max}}$")

colormap(t_axes, "turbo")
clim(t_axes, [0, 2])
% axis(t_axes, "equal")

xlabel(t_axes, "$\theta$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$r := D \tau / h^2$", Interpreter="latex", FontSize=10.5);
grid(t_axes, "on");
legend(t_axes, h, Box="off", Location="best", Interpreter="latex", FontSize=10.5)

title(t_TCL, {"(SI\_CD) $r \theta (c_{j-1}^{n+1} + c^{n+1}_{j+1}) + (1 + 2r \theta) c_j^{n+1} = r(1 - \theta) (c^n_{j-1} + c^n_{j+1}) + (1 - 2r (1 - \theta))$"; ...
              }, ...
              Interpreter="latex", FontSize=10.5);

t_txt_box = annotation(t_fig, "textbox", String="$$G = \frac{1 - 2r(1 - \theta) (1 - \cos{(\xi h)})}{1 + 2r \theta (1 - \cos{(\xi h)})}, \quad |\xi h| \le \pi$$", Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4])*0.7, .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
UNIT_ORIGINAL = t_txt_box.Units;
t_txt_box.Units = "points";
t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*0, 10.5*5], 10.5*10, 10.5*5];
t_txt_box.Units = UNIT_ORIGINAL;

if true
    print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-vector", "-dsvg")
end

%% Function:

function abs_max = G_abs_max(r, theta)
    x = 0:.001:2;
    G = (1 - 2*r*(1 - theta)*x) ./ (1 + 2*r*theta*x);
    abs_max = max(abs(G));
end

%% Function: Unit

function [f_list, x_grid, t_list] = hw4_1_unit(scheme_name, theta_list, delta_x, r_num, bndry_func_type, bndry_func, bndry_func_disp_str, init_func, init_func_disp_str, flag_save_fig, fig_name)
    arguments
        scheme_name 
        theta_list 
        delta_x
        r_num
        bndry_func_type 
        bndry_func 
        bndry_func_disp_str
        init_func
        init_func_disp_str
        flag_save_fig = false;
        fig_name = "fig_unknown"
    end

    if scheme_name == "time_SI_space_CD"
        scheme_name_short = "SI\_CD";
    end
    
    solver = DiffuProb1D();

    %%% Solver params (1-D Advection problem)

    diffu_coeff = .1;
    t_start = 0;
    x_range = [-1, 1] * (.05 / delta_x);
    % init_func = @(x) double(abs(x) < .3);
    % init_func_disp_str = "$p(x) = 1, \; |x| < 0.3$";
    % bndry_func_type = "Neumann";
    % bndry_func = @(x, t) zeros(size(x));
    % bndry_func_disp_str = "$q(x) = 0$";
    % delta_x = .05;
    % r_num = .45;
    
    delta_t = r_num*delta_x^2 / diffu_coeff;
    t_query = [0, delta_t, .1, 1, 100];  

    %%% Solving

    f_list = cell(length(theta_list), 1);
    x_grid = f_list;
    t_list = x_grid;
    for test_idx = 1:length(theta_list)
        params_struct = DiffuProb1D.prepare_params(diffu_coeff, t_start, x_range, init_func, bndry_func_type, bndry_func, delta_t, delta_x, t_query);
        solver.reset(params_struct);
        solver.solve(scheme_name, theta=theta_list(test_idx));
        [f_list{test_idx}, x_grid{test_idx}, t_list{test_idx}, ~] = solver.get_solution();
    end

    %%% Create figure

    t_fig = figure(Name=fig_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 16, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, floor(length(theta_list) / 2), 2, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$", Interpreter="latex", FontSize=10.5);
    if bndry_func_type == "Neumann"
        ylabel(t_TCL, "(numerical, " + scheme_name_short + ") $f$ (subjected to $f_{t} = Dc_{xx}, \; f|_{t = 0} = p(x), \; f_x|_{x \in \partial \Omega} = q(x)$)", Interpreter="latex", FontSize=10.5);
    elseif bndry_func_type == "Dirichlet"
        ylabel(t_TCL, "(numerical, " + scheme_name_short + ") $f$ (subjected to $f_{t} = Dc_{xx}, \; f|_{t = 0} = p(x), \; f|_{x \in \partial \Omega} = q(x)$)", Interpreter="latex", FontSize=10.5);
    else
        ylabel(t_TCL, "(numerical, " + scheme_name_short + ") $f$ (subjected to $f_{t} = Dc_{xx}, \; f|_{t = 0} = p(x), \; $ unknown B.C.)", Interpreter="latex", FontSize=10.5);
    end
    
    if scheme_name == "time_SI_space_CD"
        title(t_TCL, sprintf("%s, $\\quad$ %s, $\\quad$ $h$ = %g, $\\, \\tau$ = %g, $\\, r$ = %g", init_func_disp_str, bndry_func_disp_str, delta_x, delta_t, r_num), Interpreter="latex", FontSize=10.5);
    end

    for test_idx = 1:length(theta_list)
        t_axes = nexttile(t_TCL, test_idx);
        hold(t_axes, "on");
        for t_query_idx = 1:length(t_list{test_idx})
            plot(t_axes, x_grid{test_idx}, f_list{test_idx}{t_query_idx}, LineWidth=1.0, DisplayName=sprintf("$t = %g$", t_list{test_idx}(t_query_idx)))
        end
        set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
            Tag=sprintf("(%c) $\\theta = %g$", char('a' - 1 + test_idx), theta_list(test_idx)));
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

    if flag_save_fig
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-vector", "-dsvg")
    end
end
