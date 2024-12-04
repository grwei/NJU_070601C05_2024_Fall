%% DoSag2D.m
% Description: Solving the Streeter-Phelps model in a 2-D river
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Nov. 17, 2024
% Last modified: Dec. 4, 2024
%

%%

classdef DoSag2D < handle
    properties (Constant)
        PARAMS_REQUIRED = ["t_start", "range_x", "range_y", "velocity", "diffu_coeff", "K_r", "Q_B_loc", "Q_B", "K_a", "O_s", "init_func", "bndry_func_x0", "bndry_func_x1", "bndry_func_y0", "bndry_func_y1", "bndry_func_type_x0", "bndry_func_type_x1", "bndry_func_type_y0", "bndry_func_type_y1", "delta_t", "delta_x", "delta_y", "t_query"];
        SCHEME_NAME = ["time_forward_Euler_space_central"];
    end

    properties (Access = public)
        desc_str
    end

    properties (GetAccess = public, SetAccess = immutable)
        objID
        LOG_FILE_PATH
    end

    methods (Static, Access=public)
        params_struct = prepareParams(t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query);
    end

    methods (Access=public)
        function obj = DoSag2D(desc_str, varargin)
            arguments
                desc_str string = "description of this solver\n"
            end
            arguments (Repeating)
                varargin
            end

            obj.objID = "DoSag_" + string(datetime("now"), "yyyyMMdd_HHmmss_SS", "en_US");
            obj.desc_str = desc_str;

            if ~isfolder(".\log\")
                mkdir(".\log\")
            end
            obj.LOG_FILE_PATH = sprintf(".\\log\\%s.log", obj.objID);
            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tSolver created.\n\tDescription: %s\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), desc_str);
            fclose(LOG_FILE_ID);
            
            if nargin < 2
                obj.reset(DoSag2D.prepareParams());
                return
            end
            obj.reset(DoSag2D.prepareParams(varargin{:}));
        end

        obj = reset(obj, params_struct);
    
        [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = getSolution(obj);

        function file_path = saveSnapShot(obj, file_path)
            arguments
                obj
                file_path = ".\lib\" + obj.objID + ".mat";
            end
            [filepath, ~, ~] = fileparts(file_path);
            if ~isfolder(filepath)
                mkdir(filepath);
            end
            save(file_path, "obj");
        end

        function obj = solve(obj, scheme_name)
            arguments (Input)
                obj 
                scheme_name string = "time_forward_Euler_space_central";
            end

            assert(length(obj.velocity) == 2, "仅支持 2-D 问题");
            [obj.f_list, obj.x_grid, obj.y_grid, obj.t_list] = obj.solve2D(scheme_name);
            obj.saveSnapShot();
        end
    end

    methods (Access=private)
        f_steady = solve2DSteady(obj, scheme_name);

        [f_list, x_grid, y_grid, t_list] = solve2D(obj, scheme_name);

        B_src = bSrc(obj);
    end % end of `methods (Access=private)`

    properties (GetAccess = public, SetAccess=private)
        % solver params
        t_start
        range_x
        range_y
        velocity % [u, v]
        diffu_coeff % [D_x, D_y]
        K_r
        Q_B_loc % [x1,y1; x2, y2;] BOD 点源位置
        Q_B     % BOD 点源强度. [Q1; Q2;]
        K_a
        O_s
        init_func % {B_func; O_func}
        bndry_func_x0 % x = range_x(1) 边界 {B_func; O_func}
        bndry_func_x1 % x = range_x(2) 边界 right boundary
        bndry_func_y0 % y = range_y(1) 边界
        bndry_func_y1 % y = range_y(2) 边界
        bndry_func_type_x0 % x = range_x(1) 边界类型. "Dirichlet" | "Neumann"
        bndry_func_type_x1 % x = range_x(2) 边界类型.
        bndry_func_type_y0 % y = range_y(1) 边界类型.
        bndry_func_type_y1 % y = range_y(2) 边界类型.
        delta_t
        delta_x
        delta_y
        t_query

        % 中间结果
        B_src   % (i, j) 元素表示 x = x_grid{i}, y = y_grid{j} 格点上的 BOD 源项 (加载率): \sum_{(x_0, y_0)} Q(x_0, y_0) I_{i,j}(x_0, y_0) / (h_x h_y), 即该格点所代表的那个有限体积元中的源
        f_steady  % f_steady{1}, f_steady{2} 分别表示 BOD, DO 的稳态解

        % results
        x_grid
        y_grid
        t_list
        f_list  % f_list{k}{1}(i, j)、f_list{k}{2}(i, j) 分别是 t_list(k) 时刻 x = x_grid(i), y = y_grid(j) 格点的 BOD、DO 值.
        scheme_name
    end
end % end of `classdef DoSag2D < handle`
