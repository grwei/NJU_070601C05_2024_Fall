%% hw2_3.m
% Description: Solving the 2-D advection problem: c_t + u c_x + v c_y = 0
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 24, 2024
% Last modified: Oct. 26, 2024
%

clc; clear; close all

%% Problem 3(e)

Cr_x = 0:1e-2:2;
Cr_y = 0:1e-2:2;

[Cr_X, Cr_Y] = ndgrid(Cr_x, Cr_y);
[G_max, kh_x, kh_y] = G_max_guess(Cr_X, Cr_Y);

t_fig = figure(Name="Cr_stability");
    
% set figure size
UNIT_ORIGINAL = t_fig.Units;
t_fig.Units = "centimeters";
t_fig.Position = [3, 3, 16, 16];
t_fig.Units = UNIT_ORIGINAL;
t_TCL = tiledlayout(t_fig, 2, 2, TileSpacing="compact", Padding="tight");

% title(t_TCL, "title", Interpreter="latex", FontSize=10.5);

%%% (a) G_max

t_axes = nexttile(t_TCL, 1);

contour(t_axes, Cr_x, Cr_y, abs(G_max).', LevelList=[1,1], LineWidth=2.0)
hold on
[C, h] = contour(t_axes, Cr_x, Cr_y, abs(G_max).', LevelList=0:.25:2, LineWidth=1.0);
hold off
clabel(C, h, [.5, 1, 1.5], FontSize=8, Interpreter="latex");

set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", Tag="$|G|_{\mathrm{max}}$")
colormap(t_axes, "turbo")
clim(t_axes, [0, 2])
axis(t_axes, "equal")

xlabel(t_axes, "$\mathrm{Cr}_x$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$\mathrm{Cr}_y$", Interpreter="latex", FontSize=10.5);

% cb = colorbar(nexttile(t_TCL, 1), "eastoutside");
% cb.Layout.Tile = "east";
% cb.TickLabelInterpreter = "latex";

%%% (b) CTU example

scheme_name = "time_forward_Euler_space_CTU";
t_start = 0;
x_range = {[0, 1], [0, 1]};
init_func = @(x, y) (((x - 1/2).^2 + (y - 1/2).^2)*16 < 1) .* (1 + cos(pi * 4*sqrt( (x - 1/2).^2 + (y - 1/2).^2 ))) / 2;
bndry_func = @(x, y, t) "periodic";
velocity = [-1, 1];

% Solver params (2-D Advection problem)

delta_t = 1e-3;
Cr = [1, 1]; % Cr = u delta_t / delta_x
delta_x = abs(velocity) * delta_t ./ Cr;
t_query = 1;

solver = AdvecProb2D();
params_struct = AdvecProb2D.prepare_params(velocity, 0, x_range, init_func, bndry_func, delta_t, delta_x, t_query);
solver.reset(params_struct);
solver.solve(scheme_name);
[f_list, x_grid, t_list, ~] = solver.get_solution();

t_axes = nexttile(t_TCL, 2);

[C, h] = contour(t_axes, x_grid{1}, x_grid{2}, f_list{1}.', LevelList=0:.1:1, LineWidth=1.0);
clabel(C, h, [.1, .9], FontSize=8, Interpreter="latex");

set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex");

colormap(t_axes, "turbo")
clim(t_axes, [0, 1])
axis(t_axes, "equal")

xlabel(t_axes, "$x$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$y$", Interpreter="latex", FontSize=10.5);

%%% (c) kh_x

t_axes = nexttile(t_TCL, 3);
t_query_idx = length(t_query);

[C, h] = contour(t_axes, Cr_x, Cr_y, kh_x.', LevelList=-pi:pi/4:pi, LineWidth=1.0);
clabel(C, h, [-pi, -pi/2, 0, pi/2, pi], FontSize=8, Interpreter="latex");
set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
    Tag="$\xi_x h_x$");
colormap(t_axes, "turbo")
clim(t_axes, [-pi, pi])
axis(t_axes, "equal")

xlabel(t_axes, "$\mathrm{Cr}_x$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$\mathrm{Cr}_y$", Interpreter="latex", FontSize=10.5);

%%% (d) kh_y

t_axes = nexttile(t_TCL, 4);
t_query_idx = length(t_query);

[C, h] = contour(t_axes, Cr_x, Cr_y, kh_y.', LevelList=-pi:pi/4:pi, LineWidth=1.0);
clabel(C, h, [-pi, -pi/2, 0, pi/2, pi], FontSize=8, Interpreter="latex");
set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
    Tag="$\xi_y h_y$");
colormap(t_axes, "turbo")
clim(t_axes, [-pi, pi])
axis(t_axes, "equal")

xlabel(t_axes, "$\mathrm{Cr}_x$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$\mathrm{Cr}_y$", Interpreter="latex", FontSize=10.5);

%%% final

str_{1} = sprintf("{\\bf (b)} Cr = (%.2f, %.2f), $t$ = %.2f,", Cr(1), Cr(2), t_list);
str_{2} = sprintf(" $\\tau$ = %.1e, $\\Delta x$ = (%.1e, %.1e)", delta_t, delta_x(1), delta_x(2));
t_axes = nexttile(t_TCL, 2);
t_txt_box = annotation(t_fig, "textbox", String=str_, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
UNIT_ORIGINAL = t_txt_box.Units;
t_txt_box.Units = "points";
t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*25, 10.5*3], 10.5*25, 10.5*3];
t_txt_box.Units = UNIT_ORIGINAL;

q_ = quantile(f_list{1}, [.01, .75, .99], "all");
t_txt_box = annotation(t_fig, "textbox", String={"$q(.01,.75,.99)$ = ", sprintf("(%.1e, %.1e, %.1e)", q_(1), q_(2), q_(3))}, Position=[t_axes.Position([1, 2]) + [t_axes.Position(3), 0], .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
UNIT_ORIGINAL = t_txt_box.Units;
t_txt_box.Units = "points";
t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*25, 0], 10.5*25, 10.5*3];
t_txt_box.Units = UNIT_ORIGINAL;

for i = [1,3,4]
    t_axes = nexttile(t_TCL, i);
    t_txt_box = annotation(t_fig, "textbox", String=sprintf("{\\bf (%c)} %s", char('a'-1+i), t_axes.Tag), Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
    UNIT_ORIGINAL = t_txt_box.Units;
    t_txt_box.Units = "points";
    t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*25, 10.5*3], 10.5*25, 10.5*3];
    t_txt_box.Units = UNIT_ORIGINAL;
end

if true
    print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")
end

%% Problem 3(b)

clear;

scheme_name = "time_forward_Euler_space_upwind";
velocity = [-1, 1];

%%% initial value

delta_t_list = 1e-4;
delta_x_list = {[1e-4, 1e-4]};
t_query = 0;
flag_save = true;
fig_name_surfix = "initial";
hw2_3_unit(scheme_name, [0, 0], delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%%% fixed delta_x

delta_t_list = [.25, .05, .5, .51] * .025 / 5; % Cr = u delta_t / delta_x
delta_x_list = {[.025, .025]/5, [.025, .025]/5, [.025, .025]/5, [.025, .025]/5};
t_query = 1 * 5;
flag_save = true;
fig_name_surfix = "fixed_delta_x";
hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%%% fixed delta_t

delta_t_list = repmat(1e-3, [1, 4]);
Cr_list = {[.25, .25], [.5, .25], [.25, .5], [.5, .5]}; % Cr = u delta_t / delta_x
delta_x_list = cellfun(@(Cr, delta_t) abs(velocity) * delta_t ./ Cr, Cr_list, num2cell(delta_t_list), UniformOutput=false);
t_query = 5;
flag_save = true;
fig_name_surfix = "fixed_delta_t";
hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%%% fixed delta_t_fine

delta_t_list = repmat(2.5e-4, [1, 4]);
Cr_list = {[.25, .25], [.5, .25], [.25, .5], [.5, .5]}; % Cr = u delta_t / delta_x
delta_x_list = cellfun(@(Cr, delta_t) abs(velocity) * delta_t ./ Cr, Cr_list, num2cell(delta_t_list), UniformOutput=false);
t_query = 5;
flag_save = true;
fig_name_surfix = "fixed_delta_t_fine";
hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%% Problem 3(d)

scheme_name = "time_forward_Euler_space_CTU";
velocity = [-1, 1];

%%% fixed delta_x

Cr_list = {[.5, .5], [.25, .25], [.75, .75], [1.01, 1.01]}; % Cr = u delta_t / delta_x
delta_x_list = {[.025, .025]/5, [.025, .025]/5, [.025, .025]/5, [.025, .025]/5};
delta_t_list = cellfun(@(delta_x, Cr) Cr(1) * delta_x(1) / abs(velocity(1)), delta_x_list, Cr_list);
t_query = 5;
flag_save = true;
fig_name_surfix = "fixed_delta_x";
[f_list, x_grid, t_list] = hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%%% fixed delta_t

delta_t_list = repmat(2e-3, [1, 4]);
Cr_list = {[.5, .5], [.75, .25], [.25, .75], [1, 1]}; % Cr = u delta_t / delta_x
delta_x_list = cellfun(@(Cr, delta_t) abs(velocity) * delta_t ./ Cr, Cr_list, num2cell(delta_t_list), UniformOutput=false);
t_query = 5;
flag_save = true;
fig_name_surfix = "fixed_delta_t";
hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%%% fixed delta_t_fine

delta_t_list = repmat(5e-4, [1, 4]);
Cr_list = {[.5, .5], [.75, .25], [.25, .75], [1, 1]}; % Cr = u delta_t / delta_x
delta_x_list = cellfun(@(Cr, delta_t) abs(velocity) * delta_t ./ Cr, Cr_list, num2cell(delta_t_list), UniformOutput=false);
t_query = 5;
flag_save = true;
fig_name_surfix = "fixed_delta_t_fine";
hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix);

%% Function

function [f_list, x_grid, t_list] = hw2_3_unit(scheme_name, velocity, delta_t_list, delta_x_list, t_query, flag_save, fig_name_surfix)
    arguments (Input)
        scheme_name = "time_forward_Euler_space_upwind";
        velocity = [1, 1];
        delta_t_list = [.05, .25, .5, .51] * .025 / 5;
        delta_x_list = {[.025, .025]/5, [.025, .025]/5, [.025, .025]/5, [.025, .025]/5};
        t_query {isscalar} = 1 * 5; % must be scalar
        flag_save = false;
        fig_name_surfix = "fixed_delta_x"
    end

    solver = AdvecProb2D();

    %%% Solver params (2-D Advection problem)

    t_start = 0;
    x_range = {[0, 1], [0, 1]};
    init_func = @(x, y) (((x - 1/2).^2 + (y - 1/2).^2)*16 < 1) .* (1 + cos(pi * 4*sqrt( (x - 1/2).^2 + (y - 1/2).^2 ))) / 2;
    bndry_func = @(x, y, t) "periodic";

    %%% Solving

    % Cr_list = {[.05, .05], [.25, .25], [.5, .5], [.51, .51]}; % Cr = u delta_t / delta_x
    f_list = cell(length(delta_t_list), 1);
    x_grid = f_list;
    t_list = x_grid;
    for test_idx = 1:length(delta_t_list)
        delta_t = delta_t_list(test_idx);
        delta_x = delta_x_list{test_idx};
        params_struct = AdvecProb2D.prepare_params(velocity, t_start, x_range, init_func, bndry_func, delta_t, delta_x, t_query);
        solver.reset(params_struct);
        solver.solve(scheme_name);
        [f_list{test_idx}, x_grid{test_idx}, t_list{test_idx}, ~] = solver.get_solution();
    end

    %%% Create figure
    switch scheme_name
        case "time_forward_Euler_space_upwind"
            fig_name_prefix = "fig_3_upwind";
            scheme_name_short = "upwind";
        case "time_forward_Euler_space_CTU"
            fig_name_prefix = "fig_3_CTU";
            scheme_name_short = "CTU";
        otherwise
            fig_name_prefix = "fig3_undef";
            scheme_name_short = "unknown scheme";
    end

    t_fig = figure(Name=sprintf("%s_Cr_compare_%s", fig_name_prefix, fig_name_surfix));

    % plot
    if length(delta_t_list) > 1
        % set figure size
        UNIT_ORIGINAL = t_fig.Units;
        t_fig.Units = "centimeters";
        t_fig.Position = [3, 3, 17, 16];
        t_fig.Units = UNIT_ORIGINAL;
        t_TCL = tiledlayout(t_fig, floor(length(delta_t_list) / 2), 2, TileSpacing="compact", Padding="tight");
    else
        % set figure size
        UNIT_ORIGINAL = t_fig.Units;
        t_fig.Units = "centimeters";
        t_fig.Position = [3, 3, 12, 11];
        t_fig.Units = UNIT_ORIGINAL;
        t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="tight");
    end
    xlabel(t_TCL, "$x$", Interpreter="latex", FontSize=10.5);
    ylabel(t_TCL, "$y$", Interpreter="latex", FontSize=10.5);
    title(t_TCL, sprintf("$\\rho_t + u \\cdot \\nabla \\rho = 0$, $u$ = (%g, %g), (%s scheme)", velocity(1), velocity(2), scheme_name_short), Interpreter="latex", FontSize=10.5);
    for test_idx = 1:length(delta_t_list)
        Cr = abs(velocity) * delta_t_list(test_idx) ./ delta_x_list{test_idx};
        t_axes = nexttile(t_TCL, test_idx);
        hold(t_axes, "on");
        t_query_idx = length(t_query);
        if scheme_name == "time_forward_Euler_space_upwind" && sum(Cr) > 1 && t_query > 0 || scheme_name == "time_forward_Euler_space_CTU" && any(Cr > 1) && t_query > 0
            [C, h] = contour(t_axes, x_grid{test_idx}{1}, x_grid{test_idx}{2}, f_list{test_idx}{t_query_idx}.', 3, LineWidth=1.0, DisplayName=sprintf("$t = %.2f, \\, \\mathrm{Cr} = (%.2f, %.2f)$", t_list{test_idx}(t_query_idx), Cr(1), Cr(2)));
            clabel(C, h, FontSize=8, Interpreter="latex");
        else
            [C, h] = contour(t_axes, x_grid{test_idx}{1}, x_grid{test_idx}{2}, f_list{test_idx}{t_query_idx}.', LevelList=0:.1:1, LineWidth=1.0, DisplayName=sprintf("$t = %.2f, \\, \\mathrm{Cr} = (%.2f, %.2f)$", t_list{test_idx}(t_query_idx), Cr(1), Cr(2)));
            clabel(C, h, [.1, .9], FontSize=8, Interpreter="latex");
        end
        set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", ...
            Tag=sprintf("(%c) $\\mathrm{Cr} = (%.2f, %.2f), \\, \\tau = \\rm{%.2e}, \\, t = %.2f$", char('a' - 1 + test_idx), Cr(1), Cr(2), delta_t_list(test_idx), t_list{test_idx}(t_query_idx)));
        colormap(t_axes, "turbo")
        clim(t_axes, [0, 1])
        axis(t_axes, "equal")

        if ~mod(test_idx, 2)
            t_axes.YTickLabel = {};
        end
        if test_idx < length(delta_t_list) - 1
            t_axes.XTickLabel = {};
        end
    end

    cb = colorbar(nexttile(t_TCL, length(delta_t_list)), "east");
    cb.Layout.Tile = "east";
    cb.TickLabelInterpreter = "latex";

    % for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
    for test_idx = 1:length(delta_t_list)
        Cr = abs(velocity) * delta_t_list(test_idx) ./ delta_x_list{test_idx};
        str_{1} = sprintf("{\\bf (%c)} Cr = (%.2f, %.2f), $t$ = %.2f,", char('a'-1+test_idx), Cr(1), Cr(2), t_list{test_idx});
        str_{2} = sprintf(" $\\tau$ = %.1e, $\\Delta x$ = (%.1e, %.1e)", delta_t_list(test_idx), delta_x_list{test_idx}(1), delta_x_list{test_idx}(2));
        t_axes = nexttile(t_TCL, test_idx);
        t_txt_box = annotation(t_fig, "textbox", String=str_, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*25, 10.5*3], 10.5*25, 10.5*3];
        t_txt_box.Units = UNIT_ORIGINAL;
        if scheme_name == "time_forward_Euler_space_upwind" && sum(Cr) > 1 && t_query > 0
            t_txt_box.Color = "#4DBEEE";
        end
        if scheme_name == "time_forward_Euler_space_CTU" && any(Cr > 1) && t_query > 0
            t_txt_box.Color = "#4DBEEE";
        end

        q_ = quantile(f_list{test_idx}{1}, [.01, .75, .99], "all");
        t_txt_box = annotation(t_fig, "textbox", String={"$q(.01,.75,.99)$ = ", sprintf("(%.1e, %.1e, %.1e)", q_(1), q_(2), q_(3))}, Position=[t_axes.Position([1, 2]) + [t_axes.Position(3), 0], .1, .1], FontSize=10.5, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10.5*25, 0], 10.5*25, 10.5*3];
        t_txt_box.Units = UNIT_ORIGINAL;

        if scheme_name == "time_forward_Euler_space_upwind" && sum(Cr) > 1 && t_query > 0
            t_txt_box.Color = "#4DBEEE";
        end
        if scheme_name == "time_forward_Euler_space_CTU" && any(Cr > 1) && t_query > 0
            t_txt_box.Color = "#4DBEEE";
        end
    end

    if flag_save
        print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")
    end
end

function [G_max, kh_x_max, kh_y_max] = G_max_guess(x, y)
    arguments
        x
        y
    end
    
    size_x = size(x);
    size_y = size(y);

    if any(size_x ~= size_y)
        error("size(x) must be equal to size(y).")
    end

    k_hx = linspace(-pi,pi,100);
    k_hy = k_hx;

    [k_hx_grid, k_hy_grid] = ndgrid(k_hx, k_hy);
    exp_i_k_hx = exp(1i * k_hx_grid);
    exp_i_k_hy = exp(1i * k_hy_grid);
    [G_max, I] = arrayfun(@(x_, y_) max((1 - x_ - y_ + x_*y_) + x_*(1-y_)*exp_i_k_hx + y_*(1-x_)*exp_i_k_hy + x_*y_*exp_i_k_hx.*exp_i_k_hy, [], "all", "linear", "ComparisonMethod", "abs"), x, y);

    kh_x_max = k_hx_grid(I);
    kh_y_max = k_hy_grid(I);
end
