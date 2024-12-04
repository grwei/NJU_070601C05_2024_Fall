%% hw5_1.m
% Description: Solving the Streeter-Phelps model in a 2-D river
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Dec. 1, 2024
% Last modified: Dec. 5, 2024
%

clc; clear; close all

%%% 程序参数

FLAG_EXPORT_GRAPH = true;
RESULTS_PATH = ".\lib\hw5_1_20241203_080420.mat";
LOG_FILE_PATH = ".\log\log_" + string(datetime("now"), "yyyyMMdd_HHmmss", "en_US") + ".log";

%%% 环境准备

if ~isfolder(".\log\")
    mkdir(".\log\")
end
LOG_FILE_ID = fopen(LOG_FILE_PATH,'a');
fprintf(LOG_FILE_ID, "%s\n\tStart.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"));

%% 1. 求解

if ~isfile(RESULTS_PATH)
    t_start = tic;
    warning("%s\n\tresults file:\n\t%s\n\tdoes not exist, 正在重新求解.\n\t求解完成后, 请相应更改 script 中的 RESULTS_PATH.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH)
    fprintf(LOG_FILE_ID, "%s\n\tresults file:\n\t%s\n\tdoes not exist, 正在重新求解.\n\t求解完成后, 请相应更改 script 中的 RESULTS_PATH.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH);
    RESULTS_PATH = hw5_1_4();
    warning("%s\n\t已重新求解, 耗时 %.1f min. 请将 script 中的 RESULTS_PATH 更改为:\n\t%s\n", toc(t_start) / 60, string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH)
    fprintf(LOG_FILE_ID, "%s\n\t已重新求解, 耗时 %.1f min. 请将 script 中的 RESULTS_PATH 更改为:\n\t%s\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), toc(t_start) / 60, RESULTS_PATH);
end

%% 2. 绘图

load(RESULTS_PATH, "solvers");

time_inds = [1,2,5,6];
test_names = [
    "single-point";
    "two-points";
    "line-distributed";
    "face-distributed"
];
var_names = [
    "BOD";
    "DO";
];

%%% 2.2 各试验, 河道中心线上的 BOD, DO 稳态分布

fig_name = sprintf("fig_2");
varSteadyDistYline(solvers, 15, fig_name, FLAG_EXPORT_GRAPH);

%%% 2.1 各试验的 BOD, DO 时空演变图

for var_ind = 1:length(var_names)
    for test_ind = 1:length(solvers)
        title_str = sprintf("%s (%s discharge)", var_names(var_ind), test_names(test_ind));
        fig_name = sprintf("fig_1_%s_%s", var_names(var_ind), test_names(test_ind));
        varSpatialTemporalDist(solvers{test_ind}, var_ind, time_inds, title_str, fig_name, FLAG_EXPORT_GRAPH);
    end
end

%% 3. 生成视频

hw5_1_movie();

%% 4. 程序结束

fprintf(LOG_FILE_ID, "%s\n\tFinished.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"));
fclose(LOG_FILE_ID);

% clc; clear; close all;

%% 

function [] = varSteadyDistYline(solvers, y_query, fig_name, flag_export_graph, fig_dir)
    arguments
        solvers 
        y_query 
        fig_name 
        flag_export_graph = false;
        fig_dir = ".\fig\"
    end
    test_names = [
        "single-point discharge";
        "two-points discharge";
        "line-distributed discharge";
        "face-distributed discharge"
    ];

    t_fig = figure(Name=fig_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 18, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    t_TCL = tiledlayout(t_fig, length(solvers), 1, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
    ylabel(t_TCL, "biochemical Oxygen demand (BOD) (mg/L)", Interpreter="latex", FontSize=10, Color="#0072BD");

    Y_IND = NaN;
    for ind_ = 1:length(solvers)
        if isnan(Y_IND)
            [~, Y_IND] = min(solvers{ind_}.y_grid - y_query, [], "ComparisonMethod", "abs");
        end

        t_axes = nexttile(t_TCL, ind_);
        plot(t_axes, solvers{ind_}.x_grid, solvers{ind_}.f_list{end}{1}(:, Y_IND), LineWidth=.8);
        hold(t_axes, "on")
        set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", ...
            Tag=sprintf("(%c) %s", char('a' - 1 + ind_), test_names(ind_)));
        xline(t_axes, [10, 70], "--", Color="#A9A9A9")
        yyaxis right
        plot(t_axes, solvers{ind_}.x_grid, solvers{ind_}.f_list{end}{2}(:, Y_IND), LineWidth=.8);
        set(t_axes, YColor = "#D95319", YLim = [5.5, 8], XLimitMethod="tight");
        hold(t_axes, "off")
        yyaxis left
        t_axes.YColor = "#0072BD";
        yyaxis right
        
        if ind_ < length(solvers)
            t_axes.XTickLabel = {};
        end
    end
    title(t_TCL, sprintf("steady state ($y$ = %.3g m)", solvers{ind_}.y_grid(Y_IND)), Interpreter="latex", FontSize=10);

    linkaxes(findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}), 'xy');

    % add right ylabel
    % cb = colorbar(t_axes, Visible="off", TickLabels={});
    % cb.Position(3) = .5*cb.Position(3);
    % cb.Layout.Tile = "east";
    annotation(t_fig, "textbox", String="dissolved Oxygen (DO) (mg/L)", Color="#D95319", FitBoxToText=false, ...
        Position=[sum(t_TCL.InnerPosition([1,3])), sum(t_TCL.InnerPosition([2,4])), .89*(sum(t_TCL.InnerPosition([2,4])) - t_TCL.InnerPosition(2)), 1.2*(sum(t_TCL.OuterPosition([1,3])) - sum(t_TCL.InnerPosition([1,3])))], ...
        FontSize=10, Rotation=-90, Interpreter="latex", LineStyle="none", HorizontalAlignment="center", VerticalAlignment="top");

    % 子图编号
    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        ind_ = t_axes.Layout.Tile;
        var_data{2} = solvers{ind_}.f_list{end}{2}(:, Y_IND); % DO
        var_data{1} = solvers{ind_}.f_list{end}{1}(:, Y_IND); % BOD
        [BOD_max, x_ind_BOD] = max(var_data{1});
        [DO_min, x_ind_DO] = min(var_data{2});
        str_{3} = sprintf("min(DO) = %.3g ($x$ = %.3g)", DO_min, solvers{ind_}.x_grid(x_ind_DO));
        str_{2} = sprintf("max(BOD) = %.3g ($x$ = %.3g)", BOD_max, solvers{ind_}.x_grid(x_ind_BOD));
        str_{1} = "\bf" + t_axes.Tag;
        t_txt_box = annotation(t_fig, "textbox", String=str_{1}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*20, 10*2], 10*20, 10*2];
        t_txt_box.Units = UNIT_ORIGINAL;

        t_txt_box = annotation(t_fig, "textbox", String=str_(2:end), Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]).*[1, 0], .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*20, 0], 10*20, 10*3];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_export_graph
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        print(t_fig, fig_dir + t_fig.Name + ".svg", "-vector", "-dsvg")
        % !很不合理的现象: 下面这行语句, 导出的图片, 与从显示器看到的不同.
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".jpg", Resolution=800, BackgroundColor="none");
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".pdf", ContentType="vector", BackgroundColor="none");
        % close(t_fig)
    end
end

%%

function [] = varSpatialTemporalDist(solver, var_ind, time_inds, title_str, fig_name, flag_export_graph, fig_dir)
    arguments
        solver      DoSag2D
        var_ind     {isscalar} = 1;
        time_inds   {isvector} = [1,2,5,6];
        title_str   string = "BOD (single-point discharge)";
        fig_name    string = "fig_1_BOD_single_point";
        flag_export_graph = false;
        fig_dir     string = ".\fig\"
    end

    % var_name = "var_name";
    if var_ind == 1
        var_name = "BOD";
    elseif var_ind == 2
        var_name = "DO";
    else
        error("Unknown variable name.")
    end

    % 下载并应用 colormap

    if ~isfolder(".\inc\")
        mkdir(".\inc\")
    end
    if ~isfile(".\inc\cmocean_oxy.rgb")
        system('curl "https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/cmocean_oxy.rgb" >> ".\inc\cmocean_oxy.rgb"');
    end
    if ~isfile(".\inc\cmocean_matter.rgb")
        system('curl "https://www.ncl.ucar.edu/Document/Graphics/ColorTables/Files/cmocean_matter.rgb" >> ".\inc\cmocean_matter.rgb"');
    end
    if var_name == "DO"
        col_map = uint8(readmatrix(".\inc\cmocean_oxy.rgb", FileType="text", NumHeaderLines=2, ExpectedNumVariables=3));
        cb_str = "dissolved Oxygen (DO) (mg/L)";
        ct_levels = [6, 6.5, 7, 7.5];
        clim_ = [5.5, 8];
    elseif var_name == "BOD"
        col_map = uint8(readmatrix(".\inc\cmocean_matter.rgb", FileType="text", NumHeaderLines=2, ExpectedNumVariables=3));
        cb_str = "biochemical Oxygen demand (BOD) (mg/L)";
        ct_levels = 5;
        clim_ = [];
    else
        error("未定义 col_map");
    end

    t_fig = figure(Name=fig_name);
    colormap(t_fig, col_map)

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 18, 16];
    t_fig.Units = UNIT_ORIGINAL;

    % 决定 clim_
    if isempty(clim_)
        clim_ = [0, 1.1*max(solver.f_list{time_inds(end)}{var_ind}, [], "all")];
    end

    % plot
    t_TCL = tiledlayout(t_fig, length(time_inds), 1, TileSpacing="compact", Padding="compact");
    xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
    ylabel(t_TCL, "$y$ (m)", Interpreter="latex", FontSize=10);
    title(t_TCL, title_str, Interpreter="latex", FontSize=10);

    for ind_ = 1:length(time_inds)
        var_data = solver.f_list{time_inds(ind_)}{var_ind}.';

        t_axes = nexttile(t_TCL, ind_);
        pcolor(t_axes, solver.x_grid, solver.y_grid, var_data, FaceColor="interp", EdgeColor="none");
        hold(t_axes, "on")
        if var_name == "BOD"
            ct_levels = quantile(var_data, [.25, .50, .75, .90], "all");
        end
        [C, h] = contour(t_axes, solver.x_grid, solver.y_grid, var_data, ct_levels, LineWidth=.6, ShowText = "on", LabelFormat="%.3g", LabelColor="#808080", EdgeColor="#808080");
        clabel(C, h, FontSize=7.5, Interpreter="latex");
        hold(t_axes, "off")
        clim(t_axes, clim_);
        set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", YLimitMethod="tight", ...
            Tag=sprintf("(%c) $t$ = %.2f s", char('a' - 1 + ind_), solver.t_list(time_inds(ind_))));
        axis(t_axes, "equal")
        if ind_ < length(time_inds)
            t_axes.XTickLabel = {};
        end
    end

    cb = colorbar(t_axes, TickLabelInterpreter = "latex", FontSize = 10);
    cb.Layout.Tile = "east";
    set(cb.Label, String=cb_str, Interpreter = "latex", FontSize = 10);

    linkaxes(findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}), 'xy');
    for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
        ind_ = t_axes.Layout.Tile;
        var_data = solver.f_list{time_inds(ind_)}{var_ind};
        [min_val, min_ind_lin] = min(var_data, [], "all", "linear");
        [min_ind_x, min_ind_y] = ind2sub(size(var_data), min_ind_lin);
        [max_val, max_ind_lin] = max(var_data, [], "all", "linear");
        [max_ind_x, max_ind_y] = ind2sub(size(var_data), max_ind_lin);
        str_{1} = t_axes.Tag;
        str_{2} = sprintf("min. = %.2f at (%.2f, %.2f)\nmax. = %.2f at (%.2f, %.2f)", min_val, solver.x_grid(min_ind_x), solver.y_grid(min_ind_y), max_val, solver.x_grid(max_ind_x), solver.y_grid(max_ind_y));

        t_txt_box = annotation(t_fig, "textbox", String=str_{1}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]), .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="top");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*30, 10*1.5], 10*30, 10*1.5];
        t_txt_box.Units = UNIT_ORIGINAL;

        t_txt_box = annotation(t_fig, "textbox", String=str_{2}, Position=[t_axes.Position([1, 2]) + t_axes.Position([3, 4]).*[1, 0], .1, .1], FontSize=10, Interpreter="latex", LineStyle="none", HorizontalAlignment="right", VerticalAlignment="middle");
        UNIT_ORIGINAL = t_txt_box.Units;
        t_txt_box.Units = "points";
        t_txt_box.Position = [t_txt_box.Position([1,2]) - [10*30, 0], 10*30, 10*3];
        t_txt_box.Units = UNIT_ORIGINAL;
    end

    if flag_export_graph
        if ~isfolder(fig_dir)
            mkdir(fig_dir)
        end
        print(t_fig, fig_dir + t_fig.Name + ".svg", "-vector", "-dsvg") % 可能导致程序崩溃
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".jpg", Resolution=800, BackgroundColor="none");
        exportgraphics(t_fig, fig_dir + t_fig.Name + ".pdf", ContentType="vector", BackgroundColor="none"); % 可能导致程序崩溃
        % close(t_fig)
    end
end

%%

function results_path = hw5_1_4()
    %%% 0. 参数设置

    Q_B_SUM = 70;
    LIB_DIR = ".\lib\";

    if ~isfolder(LIB_DIR)
        mkdir(LIB_DIR);
    end

    %%% 1. 求解

    solvers = cell([1, 4]);
    solvers_path = solvers;

    % 1.1 单点排放

    desc_str = "单点排放.\n";
    Q_B_loc = [10, 15; 70, 15]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = [Q_B_SUM; 0];    % BOD 点源强度

    [solvers{1}, solvers_path{1}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.2 两点排放

    desc_str = "两点排放.\n";
    Q_B_loc = [10, 15; 70, 15]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{2}, solvers_path{2}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.3 分布排放

    desc_str = "线分布排放.\n";

    x_src = solvers{1}.x_grid(solvers{1}.x_grid >= 10 & solvers{1}.x_grid <= 70);
    Q_B_loc = [x_src, mean(solvers{1}.range_y)*ones([length(x_src), 1])]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{3}, solvers_path{3}] = createSolvers(desc_str, Q_B_loc, Q_B);

    % 1.4 面分布排放

    desc_str = "面分布排放.\n";
    x_src = solvers{1}.x_grid(solvers{1}.x_grid >= 10 & solvers{1}.x_grid <= 70);
    y_src = solvers{1}.y_grid(solvers{1}.y_grid >= 10 & solvers{1}.y_grid <= 20);
    [X, Y] = ndgrid(x_src, y_src);
    Q_B_loc = [X(:), Y(:)]; % [x1,y1; x2, y2;] BOD 点源位置
    Q_B = Q_B_SUM / size(Q_B_loc, 1) * ones([size(Q_B_loc, 1), 1]);    % BOD 点源强度

    [solvers{4}, solvers_path{4}] = createSolvers(desc_str, Q_B_loc, Q_B);

    %%% 2. 保存结果

    results_path = LIB_DIR + sprintf("hw5_1_%s.mat", string(datetime("now"), "yyyyMMdd_HHmmss", "en_US"));
    save(results_path, "solvers");
    delete(solvers_path{:});
end

%% 求解

function [solver, solver_path] = createSolvers(desc_str, Q_B_loc, Q_B, t_query, delta_x, delta_y, delta_t, t_start, range_x, range_y, velocity, diffu_coeff, K_r, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1)
    arguments (Input)
        desc_str string = "description of this solver\n"
        Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
        Q_B = [70; 0];    % BOD 点源强度
        t_query = [0, 15.64, 41.576, 84.26, 112.036, +Inf]; % [0, 15.64, 41.576, 84.26, 112.036, +Inf];
        delta_x = .05; % .05 | .1 | .2 | 1
        delta_y = delta_x;
        delta_t = .00124; % .001249 | .00499 | .0199 | .4950
        t_start = 0;
        range_x = [0, 150];
        range_y = [0, 30];
        velocity = [.4, 0];
        diffu_coeff = [.5, .5];
        K_r = .01;
        K_a = .02;
        O_s = 8;
        init_func = {@(x, y) zeros(size(x)); @(x, y) O_s*ones(size(x))};
        bndry_func_type_x0 = "Dirichlet";
        bndry_func_x0 = {@(y, t) zeros(size(y)); @(y, t) O_s*ones(size(y))};
        bndry_func_type_x1 = "Neumann";
        bndry_func_x1 = {@(y, t) zeros(size(y)); @(y, t) zeros(size(y))};
        bndry_func_type_y0 = "Neumann";
        bndry_func_y0 = {@(x, t) zeros(size(x)); @(x, t) zeros(size(x))};
        bndry_func_type_y1 = bndry_func_type_y0;
        bndry_func_y1 = bndry_func_y0;
    end
    arguments (Output)
        solver DoSag2D
        solver_path {isfile}
    end

    %%% 1. 输入合法性检查

    % 1.1 FTCS 格式稳定性

    % Cr = delta_t * velocity ./ [delta_x, delta_y];
    r = diffu_coeff .* delta_t ./ [delta_x, delta_y].^2;
    if any(abs(1 - 2*sum(r) - [K_r, K_a]*delta_t) + 2*sum(r) > 1)
        warning("格式的稳定性存疑. 建议 delta_t < %g", 1 / (4*diffu_coeff(1)/delta_x^2 + max([K_a, K_r])));
    end

    %%% 2. 求解

    solver = DoSag2D(desc_str, t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query);
    solver.solve("time_forward_Euler_space_central");
    % [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = solver.getSolution();
    solver_path = solver.saveSnapShot();
end
