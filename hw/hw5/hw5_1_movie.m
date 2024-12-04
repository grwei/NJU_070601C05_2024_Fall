%% hw5_1_movie.m
% Description: Solving the Streeter-Phelps model in a 2-D river
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Dec. 4, 2024
% Last modified: Dec. 5, 2024
%

function [] = hw5_1_movie()

%%% 程序参数

RESULTS_PATH = ".\lib\hw5_1_20241204_230422.mat";
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

    str_ = sprintf("%s\n\tresults file:\n\t%s\n\tdoes not exist, 正在重新求解.\n\t求解完成后, 请相应更改 script 中的 RESULTS_PATH.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), RESULTS_PATH);
    warning(str_)
    fprintf(LOG_FILE_ID, str_);

    RESULTS_PATH = hw5_1_4();

    str_ = fprintf("%s\n\t已重新求解, 耗时 %.1f min. 请将 script 中的 RESULTS_PATH 更改为:\n\t%s\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), toc(t_start) / 60, RESULTS_PATH);
    warning(str_)
    fprintf(LOG_FILE_ID, str_);
end

%% 2. 创建视频

load(RESULTS_PATH, "solvers");

MOVIE_DIR = ".\fig\movie\";
TIME_INDS = [];
VAR_NAMES = [
    "BOD";
    "DO";
];

%%% 2.1 各试验的 BOD, DO 时空演变图

for var_ind = 1:length(VAR_NAMES)
    movie_name = sprintf("mov_1_%s", VAR_NAMES(var_ind));
    varSpatialTemporalDist(solvers, var_ind, TIME_INDS, movie_name, MOVIE_DIR);
end

%%% 2.2 各试验, 河道中心线上的 BOD, DO 分布

movie_name = sprintf("mov_2_central_line");
varSteadyDistYline(solvers, 15, TIME_INDS, movie_name, MOVIE_DIR);

%% 3. 程序结束

fprintf(LOG_FILE_ID, "%s\n\tFinished.\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"));
fclose(LOG_FILE_ID);

end % end of `function [] = hw5_1_movie()`

%% 

function [] = varSteadyDistYline(solvers, y_query, time_inds, movie_name, movie_dir)
    arguments
        solvers 
        y_query 
        time_inds = [];
        movie_name = "mov_2_central_line"
        movie_dir = ".\fig\movie\"
    end
    test_names = [
        "single-point discharge";
        "two-points discharge";
        "line-distributed discharge";
        "face-distributed discharge"
    ];

    if ~isfolder(movie_dir)
        mkdir(movie_dir)
    end

    if isempty(time_inds)
        time_inds = 1:length(solvers{1}.f_list);
    end

    Y_IND = NaN;
    if isnan(Y_IND)
        [~, Y_IND] = min(solvers{1}.y_grid - y_query, [], "ComparisonMethod", "abs");
    end
    BOD_ylim_ = cellfun(@(solver) [0, 1.1*max(solver.f_list{time_inds(end)}{1}, [], "all")], solvers, "UniformOutput", false);

    h_vw = VideoWriter(movie_dir + movie_name, 'Motion JPEG AVI');
    h_vw.Quality = 100;
    h_vw.FrameRate = 10;
    h_vw.open();

    t_fig = figure(Name=movie_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 18, 16];
    t_fig.Units = UNIT_ORIGINAL;
    for ind_ = 1:length(time_inds)
        % plot
        t_TCL = tiledlayout(t_fig, length(solvers), 1, TileSpacing="compact", Padding="compact");
        title(t_TCL, sprintf("$t$ = %.2f s ($y$ = %.3g m)", solvers{1}.t_list(time_inds(ind_)), solvers{1}.y_grid(Y_IND)), Interpreter="latex", FontSize=10);
        xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
        ylabel(t_TCL, "biochemical Oxygen demand (BOD) (mg/L)", Interpreter="latex", FontSize=10, Color="#0072BD");

        for sol_ind_ = 1:length(solvers)
            t_axes = nexttile(t_TCL, sol_ind_);
            plot(t_axes, solvers{sol_ind_}.x_grid, solvers{sol_ind_}.f_list{time_inds(ind_)}{1}(:, Y_IND), LineWidth=.8);
            hold(t_axes, "on")
            set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", YLim=BOD_ylim_{sol_ind_}, ...
                Tag=sprintf("(%c) %s", char('a' - 1 + sol_ind_), test_names(sol_ind_)));
            xline(t_axes, [10, 70], "--", Color="#A9A9A9")
            yyaxis right
            plot(t_axes, solvers{sol_ind_}.x_grid, solvers{sol_ind_}.f_list{time_inds(ind_)}{2}(:, Y_IND), LineWidth=.8);
            set(t_axes, YColor = "#D95319", YLim = [5.5, 8], XLimitMethod="tight");
            hold(t_axes, "off")
            yyaxis left
            t_axes.YColor = "#0072BD";
            yyaxis right
            
            if sol_ind_ < length(solvers)
                t_axes.XTickLabel = {};
            end
        end
    
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
            sol_ind_ = t_axes.Layout.Tile;
            var_data{2} = solvers{sol_ind_}.f_list{time_inds(ind_)}{2}(:, Y_IND); % DO
            var_data{1} = solvers{sol_ind_}.f_list{time_inds(ind_)}{1}(:, Y_IND); % BOD
            [BOD_max, x_ind_BOD] = max(var_data{1});
            [DO_min, x_ind_DO] = min(var_data{2});
            str_{3} = sprintf("min(DO) = %.3g ($x$ = %.3g)", DO_min, solvers{sol_ind_}.x_grid(x_ind_DO));
            str_{2} = sprintf("max(BOD) = %.3g ($x$ = %.3g)", BOD_max, solvers{sol_ind_}.x_grid(x_ind_BOD));
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
        cdata = print(t_fig, '-RGBImage','-r600');
        F = im2frame(cdata);
        h_vw.writeVideo(F);
        clf(t_fig, "reset");
    end % end of `for ind_ = 1:length(time_inds)`
    h_vw.close();
    close(t_fig);
end

%%

function [] = varSpatialTemporalDist(solvers, var_ind, time_inds, movie_name, movie_dir)
    arguments
        solvers     
        var_ind     {isscalar} = 1;
        time_inds   = [];
        movie_name    string = "mov_1_BOD";
        movie_dir     string = ".\fig\movie\"
    end

    TEST_NAMES = [
        "single-point discharge";
        "two-points discharge";
        "line-distributed discharge";
        "face-distributed discharge"
    ];

    if ~isfolder(movie_dir)
        mkdir(movie_dir)
    end

    % var_name = "var_name";
    if var_ind == 1
        var_name = "BOD";
        % title_str = "biochemical Oxygen demand (BOD)";
    elseif var_ind == 2
        var_name = "DO";
        % title_str = "dissolved Oxygen (DO) (mg/L)";
    else
        error("Unknown variable name.")
    end

    if isempty(time_inds)
        time_inds = 1:length(solvers{1}.f_list);
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
        clim_ = {[5.5, 8]};
    elseif var_name == "BOD"
        col_map = uint8(readmatrix(".\inc\cmocean_matter.rgb", FileType="text", NumHeaderLines=2, ExpectedNumVariables=3));
        cb_str = "biochemical Oxygen demand (BOD) (mg/L)";
        ct_levels = 5;
        clim_ = cellfun(@(solver) [0, 1.1*max(solver.f_list{time_inds(end)}{var_ind}, [], "all")], solvers, "UniformOutput", false);
    else
        error("未定义 col_map");
    end

    % h_vw = VideoWriter(movie_dir + movie_name, 'Motion JPEG AVI');
    h_vw = VideoWriter(movie_dir + movie_name, 'Motion JPEG AVI');
    h_vw.Quality = 100;
    h_vw.FrameRate = 10;
    h_vw.open();

    t_fig = figure(Name=movie_name);

    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    if var_name == "DO"
        t_fig.Position = [3, 3, 18, 16];
    else
        t_fig.Position = [3, 3, 18, 16.0];
    end
    t_fig.Units = UNIT_ORIGINAL;

    % plot
    for ind_ = 1:length(time_inds)
        colormap(t_fig, col_map);
        t_TCL = tiledlayout(t_fig, length(solvers), 1, TileSpacing="compact", Padding="compact");
        xlabel(t_TCL, "$x$ (m)", Interpreter="latex", FontSize=10);
        ylabel(t_TCL, "$y$ (m)", Interpreter="latex", FontSize=10);
        title(t_TCL, sprintf("$t$ = %.2f s", solvers{1}.t_list(time_inds(ind_))), Interpreter="latex", FontSize=10);
        for sol_ind = 1:length(solvers)
            var_data = solvers{sol_ind}.f_list{time_inds(ind_)}{var_ind}.';

            t_axes = nexttile(t_TCL, sol_ind);
            pcolor(t_axes, solvers{sol_ind}.x_grid, solvers{sol_ind}.y_grid, var_data, FaceColor="interp", EdgeColor="none");
            hold(t_axes, "on")
            if var_name == "BOD"
                ct_levels = quantile(var_data, [.25, .50, .75, .90], "all");
            end
            [C, h] = contour(t_axes, solvers{sol_ind}.x_grid, solvers{sol_ind}.y_grid, var_data, ct_levels, LineWidth=.6, ShowText = "on", LabelFormat="%.3g", LabelColor="#808080", EdgeColor="#808080");
            clabel(C, h, FontSize=7.5, Interpreter="latex");
            hold(t_axes, "off")
            if isscalar(clim_)
                clim(t_axes, clim_{1});
            else
                clim(t_axes, clim_{sol_ind});
            end
            set(t_axes, FontName="Times New Roman", FontSize=10, Layer="top", Box="on", TickLabelInterpreter="latex", XLimitMethod="tight", YLimitMethod="tight", ...
                Tag=sprintf("(%c) %s", char('a' - 1 + sol_ind), TEST_NAMES(sol_ind)));
            axis(t_axes, "equal")
            if sol_ind < length(length(solvers))
                t_axes.XTickLabel = {};
            end

            if var_name ~= "DO"
                colorbar(t_axes, Location="eastoutside", TickLabelInterpreter = "latex", FontSize = 10);
            end
        end % end of `for sol_ind = 1:length(solvers)`
        if var_name == "DO"
            cb = colorbar(t_axes, Location="eastoutside", TickLabelInterpreter = "latex", FontSize = 10);
            set(cb.Label, String=cb_str, Interpreter = "latex", FontSize = 10);
            cb.Layout.Tile = "east";
        end
        linkaxes(findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}), 'xy');
        for t_axes = findobj(t_TCL, 'Type', "Axes", {'-regexp', 'Tag', "^\([a-z]+\)"}).'
            sol_ind = t_axes.Layout.Tile;
            var_data = solvers{sol_ind}.f_list{time_inds(ind_)}{var_ind};
            [min_val, min_ind_lin] = min(var_data, [], "all", "linear");
            [min_ind_x, min_ind_y] = ind2sub(size(var_data), min_ind_lin);
            [max_val, max_ind_lin] = max(var_data, [], "all", "linear");
            [max_ind_x, max_ind_y] = ind2sub(size(var_data), max_ind_lin);
            str_{1} = "\bf" + t_axes.Tag;
            str_{2} = sprintf("min. = %.2f at (%.2f, %.2f)\nmax. = %.2f at (%.2f, %.2f)", min_val, solvers{sol_ind}.x_grid(min_ind_x), solvers{sol_ind}.y_grid(min_ind_y), max_val, solvers{sol_ind}.x_grid(max_ind_x), solvers{sol_ind}.y_grid(max_ind_y));

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

        if var_name == "BOD"
            % add right ylabel
            % cb = colorbar(t_axes, Visible="off", TickLabels={});
            % cb.Position(3) = .5*cb.Position(3);
            % cb.Layout.Tile = "east";
            annotation(t_fig, "textbox", String=cb_str, Color="k", FitBoxToText=false, ...
                Position=[sum(t_TCL.InnerPosition([1,3]))*.99, sum(t_TCL.InnerPosition([2,4])), .89*(sum(t_TCL.InnerPosition([2,4])) - t_TCL.InnerPosition(2)), 1.2*(sum(t_TCL.OuterPosition([1,3])) - sum(t_TCL.InnerPosition([1,3])))], ...
                FontSize=10, Rotation=-90, Interpreter="latex", LineStyle="none", HorizontalAlignment="center", VerticalAlignment="top");
        end

        cdata = print(t_fig, '-RGBImage','-r600');
        F = im2frame(cdata);
        h_vw.writeVideo(F);
        clf(t_fig, "reset");
    end % end of `for ind_ = 1:length(time_inds)`
    h_vw.close();
    close(t_fig);
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
        t_query = [0:1:99, 100:2:200, +Inf]; % [0, 15.64, 41.576, 84.26, 112.036, +Inf];
        delta_x = .2; % .05 | .1 | .2 | 1
        delta_y = delta_x;
        delta_t = .01; % .001249 | .00499 | .0199 | .4950
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
