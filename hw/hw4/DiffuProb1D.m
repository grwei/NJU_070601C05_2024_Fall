%% DiffuProb1D.m
% Description: Solving the 1-D diffusion problem: c_t = D c_xx, D constant
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Nov. 8, 2024
% Last modified: Nov. 9, 2024
%

%%

classdef DiffuProb1D < handle
    properties (Constant)
        PARAMS_REQUIRED = ["diffu_coeff","t_start","x_range","init_func", "bndry_func_type", "bndry_func", "delta_t", "delta_x", "t_query"];
        SCHEME_NAME = ["time_SI_space_CD"];
    end

    properties (Access=private)
        % solver params
        diffu_coeff
        t_start
        x_range
        init_func
        bndry_func_type
        bndry_func
        delta_t
        delta_x
        t_query

        % results
        x_grid
        t_list
        f_list
        scheme_name

        % 三对角方程求解器
        tri_diag_solver
    end

    methods (Static, Access=public)
        function params = prepare_params(diffu_coeff, t_start, x_range, init_func, bndry_func_type, bndry_func, delta_t, delta_x, t_query)
            arguments (Input)
                diffu_coeff = .1;
                t_start = 0;
                x_range = [-1, 1];
                init_func = @(x) double(abs(x) < .3);
                bndry_func_type = "Neumann";
                bndry_func = @(x, t) 0;
                delta_t = .01125;
                delta_x = .05;
                t_query = [0, .1, 1];  
            end

            params.diffu_coeff = diffu_coeff;
            params.t_start = t_start;
            params.x_range = sort(x_range, "ascend");
            if ~iscolumn(init_func)
                init_func = init_func.';
            end
            params.init_func = init_func;
            params.bndry_func_type = bndry_func_type;
            params.bndry_func = bndry_func;
            params.delta_t = delta_t;
            params.delta_x = delta_x;
            params.t_query = sort(t_query(t_query >= t_start), "ascend");
        end
    end

    methods (Access=public)
        function obj = DiffuProb1D(varargin)
            obj.tri_diag_solver = TriDiagSolver();
            if nargin < 1
                obj.reset(DiffuProb1D.prepare_params());
                return
            end
            obj.reset(varargin{:});
        end

        function obj = reset(obj, params)
            field_missing_flag = ~isfield(params, DiffuProb1D.PARAMS_REQUIRED);
            if any(field_missing_flag)
                warning("Invalid params_struct. Params not changed. Missing field: %s", join(DiffuProb1D.PARAMS_REQUIRED(field_missing_flag), ", "));
                return
            end

            params = DiffuProb1D.prepare_params(params.diffu_coeff, params.t_start, params.x_range, params.init_func, params.bndry_func_type, params.bndry_func, params.delta_t, params.delta_x, params.t_query);
            for field_name = DiffuProb1D.PARAMS_REQUIRED
                obj.(field_name) = params.(field_name);
            end
            obj.tri_diag_solver.reset();

            obj.t_list = nan(size(obj.t_query));
            obj.f_list = cell(size(obj.t_query));
        end
    
        function [f_list, x_grid, t_list, scheme_name] = get_solution(obj)
            f_list = obj.f_list;
            x_grid = obj.x_grid;
            t_list = obj.t_list;
            scheme_name = obj.scheme_name;
        end

        function obj = solve(obj, scheme_name, opts)
            arguments (Input)
                obj 
                scheme_name
                opts.theta = .5;
            end

            if isscalar(obj.diffu_coeff)
                [obj.f_list, obj.x_grid, obj.t_list] = obj.solve_1D(scheme_name, theta=opts.theta);
                return
            end

            error("diffu_coeff must be scalar.")
        end
    end

    methods (Access=private)
        function [f_list, x_grid, t_list] = solve_1D(obj, scheme_name, opts)
            arguments (Input)
                obj
                scheme_name;
                opts.theta = .5;
            end

            f_list = cell(size(obj.t_query));
            t_list = nan(size(obj.t_query));
            x_grid = (obj.x_range(1): obj.delta_x: obj.x_range(2)).';

            %%% scheme I (2-level)

            if ismember(scheme_name, ["time_SI_space_CD"])
                obj.scheme_name = scheme_name;
                r_num = obj.diffu_coeff * obj.delta_t / obj.delta_x^2;

                f_last = obj.init_func(x_grid);
                t_last = obj.t_start;
                t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

                % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
                if isempty(t_query_next_idx) || t_query_next_idx > 1
                    t_list(1) = t_last;
                    f_list{1} = f_last;
                end
                
                % 用于记录进度
                REC_TIME_INT = 5; % seconds
                time_rec_start = tic;
                t_last_rec = t_last;

                while t_last < obj.t_query(end)
                    t_now = t_last + obj.delta_t;

                    % 输出进度
                    if toc(time_rec_start) > REC_TIME_INT
                        time_remained = REC_TIME_INT * (obj.t_query(end) - t_now) / (t_now - t_last_rec);
                        fprintf("%s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), r_num = (%s),\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), scheme_name, time_remained, join(string(size(X)), ", "), join(string(r_num), ", "), t_now, obj.t_query(t_query_next_idx), obj.t_query(end));
                        time_rec_start = tic;
                        t_last_rec = t_now;
                    end

                    % add halo points, advance one time step
                    if scheme_name == "time_SI_space_CD"
                        if obj.bndry_func_type == "periodic"
                            error("周期边界的求解程序尚未完工")
                        end
                        if obj.bndry_func_type == "Dirichlet"
                            tri_diag_coeff_mat = [[-r_num*opts.theta*ones(length(x_grid)-1, 1); 0], [1; (1 + 2*r_num*opts.theta)*ones(length(x_grid)-2, 1); 1], [0; -r_num*opts.theta*ones(length(x_grid)-1, 1)]];
                            const_vec = [obj.bndry_func(x_grid(1), t_now);
                                         r_num*(1 - opts.theta)*(f_last(1:end-2) + f_last(3:end)) + (1 - 2*r_num*(1 - opts.theta))*f_last(2:end-1);
                                         obj.bndry_func(x_grid(end), t_now)];
                            f_now = obj.tri_diag_solver.reset(tri_diag_coeff_mat, const_vec).solve();
                        end
                        if obj.bndry_func_type == "Neumann"
                            tri_diag_coeff_mat = [[-r_num*opts.theta*ones(length(x_grid)-1, 1); -2*r_num*opts.theta], (1 + 2*r_num*opts.theta)*ones(length(x_grid), 1), [-2*r_num*opts.theta; -r_num*opts.theta*ones(length(x_grid)-1, 1)]];
                            const_vec = [r_num*(1 - opts.theta)*2*(f_last(2) - obj.delta_x*obj.bndry_func(x_grid(1), t_last)) + f_last(1)*(1 - 2*r_num*(1 - opts.theta)) - 2*r_num*opts.theta*obj.bndry_func(x_grid(1), t_now);
                                         r_num*(1 - opts.theta)*(f_last(1:end-2) + f_last(3:end)) + (1 - 2*r_num*(1 - opts.theta))*f_last(2:end-1);
                                         r_num*(1 - opts.theta)*2*(f_last(end-1) + obj.delta_x*obj.bndry_func(x_grid(end), t_last)) + f_last(end)*(1 - 2*r_num*(1 - opts.theta)) + 2*r_num*opts.theta*obj.bndry_func(x_grid(end), t_now)];
                            f_now = obj.tri_diag_solver.reset(tri_diag_coeff_mat, const_vec).solve();
                        end
                    else
                        error("Unknown scheme.")
                    end

                    % store solution at queried time.
                    if ~((obj.t_query(t_query_next_idx) - t_last)*(obj.t_query(t_query_next_idx) - t_now) > 0)
                        if abs(obj.t_query(t_query_next_idx) - t_now) < abs(obj.t_query(t_query_next_idx) - t_last)
                            t_list(t_query_next_idx) = t_now;
                            f_list{t_query_next_idx} = f_now;
                        else
                            t_list(t_query_next_idx) = t_last;
                            f_list{t_query_next_idx} = f_last;
                        end
                        t_query_next_idx = t_query_next_idx + 1;
                    end

                    % prepare for next time step
                    t_last = t_now;
                    f_last = f_now;
                end % end of `while t_last < obj.t_query(end)`


                t_list_disp_str = string(t_list);
                t_list_disp_str(ismissing(t_list_disp_str)) = "NaN";
                fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(t_list_disp_str, ", "));
                return
            end % end of `if ismember(scheme_name, ["time_forward_Euler_space_upwind", "time_forward_Euler_space_central_diff"])`

            error("Invalid scheme_name. Please choose one of the following options: %s", join(DiffuProb1D.SCHEME_NAME, ", "));
        end % end of `function [f_list, x_grid, t_list] = solve_1D(obj, scheme_name)`
    end % end of `methods (Access=private)`
end % end of `classdef DiffuProb1D < handle`
