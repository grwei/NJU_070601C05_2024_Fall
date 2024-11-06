%% AdvecProb2D.m
% Description: Solving the 1-D or 2-D advection problem: c_t + u c_x + v c_y = 0, u, v constants
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 21, 2024
% Last modified: Nov. 7, 2024
%

%%

classdef AdvecProb2D < handle
    properties (Constant)
        PARAMS_REQUIRED = ["velocity","t_start","x_range","init_func","bndry_func", "delta_t", "delta_x", "t_query"];
        SCHEME_NAME = ["time_forward_Euler_space_upwind", "time_forward_Euler_space_central_diff", "time_forward_Euler_space_CTU"];
    end

    properties (Access=private)
        % solver params
        velocity
        t_start
        x_range
        init_func
        bndry_func
        delta_t
        delta_x
        t_query

        % results
        x_grid
        t_list
        f_list
        scheme_name
    end

    methods (Static, Access=public)
        function params_struct = prepare_params(velocity, t_start, x_range, init_func, bndry_func, delta_t, delta_x, t_query)
            arguments (Input)
                velocity = [1, 1];
                t_start = 0;
                x_range = {[0, 1], [0, 1]};
                init_func = @(x, y) (((x - 1/2).^2 + (y - 1/2).^2)*16 < 1) .* (1 + cos(pi * 4*sqrt( (x - 1/2).^2 + (y - 1/2).^2 ))) / 2;
                bndry_func = @(x, y, t) "periodic";
                delta_t = .0125 / 5;
                delta_x = [.0250, .0250] / 5;
                t_query = [0, .5, 1] * 5;  
            end

            params_struct.velocity = velocity;
            params_struct.t_start = t_start;
            params_struct.x_range = cellfun(@(list_) sort(list_, "ascend"), x_range, UniformOutput=false); 
            params_struct.init_func = init_func;
            params_struct.bndry_func = bndry_func;
            params_struct.delta_t = delta_t;
            params_struct.delta_x = delta_x;
            params_struct.t_query = sort(t_query(t_query >= t_start), "ascend");
        end
    end

    methods (Access=public)
        function obj = AdvecProb2D(varargin)
            if nargin < 1
                obj.reset(AdvecProb2D.prepare_params());
                return
            end
            obj.reset(varargin{:});
        end

        function obj = reset(obj, params_struct)
            field_missing_flag = ~isfield(params_struct, AdvecProb2D.PARAMS_REQUIRED);
            if any(field_missing_flag)
                warning("Invalid params_struct. Params not changed. Missing field: %s", join(AdvecProb2D.PARAMS_REQUIRED(field_missing_flag), ", "));
                return
            end

            params_struct = AdvecProb2D.prepare_params(params_struct.velocity, params_struct.t_start, params_struct.x_range, params_struct.init_func, params_struct.bndry_func, params_struct.delta_t, params_struct.delta_x, params_struct.t_query);
            for field_name = AdvecProb2D.PARAMS_REQUIRED
                obj.(field_name) = params_struct.(field_name);
            end

            obj.t_list = nan(size(obj.t_query));
            obj.f_list = cell(size(obj.t_query));
        end
    
        function [f_list, x_grid, t_list, scheme_name] = get_solution(obj)
            f_list = obj.f_list;
            x_grid = obj.x_grid;
            t_list = obj.t_list;
            scheme_name = obj.scheme_name;
        end

        function obj = solve(obj, scheme_name, ra_opts)
            arguments (Input)
                obj 
                scheme_name
                ra_opts.nu = 0;
                ra_opts.alpha = 1;
            end

            if isscalar(obj.velocity)
                [obj.f_list, obj.x_grid, obj.t_list] = obj.solve_1D(scheme_name, nu=ra_opts.nu, alpha=ra_opts.alpha);
                return
            end
            if length(obj.velocity) == 2
                [obj.f_list, obj.x_grid, obj.t_list] = obj.solve_2D(scheme_name);
                return
            end
        end
    end

    methods (Access=private)
        function [f_list, x_grid, t_list] = solve_1D(obj, scheme_name, ra_opts)
            arguments (Input)
                obj
                scheme_name;
                ra_opts.nu = 0; % Robert-Asselin filter params
                ra_opts.alpha = 1;
            end

            f_list = cell(size(obj.t_query));
            t_list = nan(size(obj.t_query));
            x_grid = cellfun(@(range_1D, delta_1D) (range_1D(1): delta_1D: range_1D(2)).', obj.x_range, num2cell(obj.delta_x), UniformOutput=false);

            %%% scheme I (2-level)

            if ismember(scheme_name, ["time_forward_Euler_space_upwind", "time_forward_Euler_space_central_diff"])
                obj.scheme_name = scheme_name;

                Cr = abs(obj.velocity) * obj.delta_t ./ obj.delta_x;
                f_last = obj.init_func(x_grid{1});
                t_last = obj.t_start;
                t_now = t_last;
                t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

                % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
                if isempty(t_query_next_idx)
                    t_list(1) = t_last;
                    f_list{1} = f_last;

                    fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                    return
                end
                
                % 用于记录进度
                REC_TIME_INT = 5; % seconds
                time_start = tic;
                t_now_last_rec = t_last;

                while t_now < obj.t_query(end)
                    t_now = t_last + obj.delta_t;

                    % 输出进度
                    
                    if toc(time_start) > REC_TIME_INT
                        time_remained = REC_TIME_INT * (obj.t_query(end) - t_now) / (t_now - t_now_last_rec);
                        fprintf("%s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), Cr = (%s),\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), scheme_name, time_remained, join(string(size(X)), ", "), join(string(Cr), ", "), t_now, obj.t_query(t_query_next_idx), obj.t_query(end));
                        time_start = tic;
                        t_now_last_rec = t_now;
                    end

                    % add halo points, advance one time step
                    if scheme_name == "time_forward_Euler_space_upwind" && obj.velocity(1) > 0
                        if string(obj.bndry_func(0, 0)) == "periodic"
                            halo_left = f_last(end-1);
                        else
                            halo_left = obj.bndry_func(x_grid{1}(1) - obj.delta_x(1), t_last);
                        end
                        halo_right = [];
                        f_last_expanded = [halo_left; f_last; halo_right];
                        f_now = (1 - Cr(1))*f_last_expanded(1 + length(halo_left) : end-length(halo_right)) + Cr(1)*f_last_expanded(length(halo_left) : -1 + end - length(halo_right));
                    elseif scheme_name == "time_forward_Euler_space_upwind"
                        halo_left = [];
                        if string(obj.bndry_func(0, 0)) == "periodic"
                            halo_right = f_last(2);
                        else
                            halo_right = obj.bndry_func(x_grid{1}(end) + obj.delta_x(1), t_last);
                        end
                        f_last_expanded = [halo_left; f_last; halo_right];
                        f_now = (1 - Cr(1))*f_last_expanded(1 + length(halo_left) : end-length(halo_right)) + Cr(1)*f_last_expanded(2 + length(halo_left) : 1 + end-length(halo_right));
                    elseif scheme_name == "time_forward_Euler_space_central_diff"
                        if string(obj.bndry_func(0, 0)) == "periodic"
                            halo_left = f_last(end-1);
                            halo_right = f_last(2);
                        else
                            halo_left = obj.bndry_func(x_grid{1}(1) - obj.delta_x(1), t_last);
                            halo_right = obj.bndry_func(x_grid{1}(end) + obj.delta_x(1), t_last);
                        end
                        f_last_expanded = [halo_left; f_last; halo_right];
                        f_now = f_last_expanded(1 + length(halo_left) : end-length(halo_right)) + sign(obj.velocity)*Cr(1)/2*(f_last_expanded(length(halo_left) : -1 + end-length(halo_right)) - f_last_expanded(2 + length(halo_left) : 1 + end-length(halo_right)));
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
                end % end of `while t_now < obj.t_query(end)`

                fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                return
            end % end of `if ismember(scheme_name, ["time_forward_Euler_space_upwind", "time_forward_Euler_space_central_diff"])`

            %%% scheme II (3-level)

            if ismember(scheme_name, ["time_leap_frog_RA_space_central_diff"])
                obj.scheme_name = scheme_name;
                Cr = obj.velocity * obj.delta_t ./ obj.delta_x;

                f_llast = obj.init_func(x_grid{1});
                t_llast = obj.t_start;
                t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

                % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
                if isempty(t_query_next_idx)
                    t_list(1) = t_llast;
                    f_list{1} = f_llast;

                    fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                    return
                end

                if t_query_next_idx > 1
                    t_list(t_query_next_idx - 1) = t_llast;
                    f_list{t_query_next_idx - 1} = f_llast;
                end

                % 先用两层方案 (upwind) 积分一个时间步
                t_last = t_llast + obj.delta_t;
                if string(obj.bndry_func(0, 0)) == "periodic"
                    halo_left = f_llast(end-1);
                    halo_right = f_llast(2);
                else
                    halo_left = obj.bndry_func(x_grid{1}(1) - obj.delta_x(1), t_llast);
                    halo_right = obj.bndry_func(x_grid{1}(end) + obj.delta_x(1), t_llast);
                end
                f_llast_expanded = [halo_left; f_llast; halo_right];
                f_last = (1 - abs(Cr(1)))*f_llast_expanded(1 + length(halo_left) : end-length(halo_right)) ...
                        + abs(Cr(1))*(obj.velocity > 0)*f_llast_expanded(length(halo_left) : -1 + end - length(halo_right)) ...
                        + abs(Cr(1))*(obj.velocity < 0)*f_llast_expanded(2 + length(halo_left) : 1 + end - length(halo_right));
                
                % 用于记录进度
                REC_TIME_INT = 5; % seconds
                time_start = tic;
                t_now_last_rec = t_last;

                % 开始积分
                while t_last < obj.t_query(end) + obj.delta_t
                    t_now = t_last + obj.delta_t;

                    % 输出进度
                    
                    if toc(time_start) > REC_TIME_INT
                        time_remained = REC_TIME_INT * (obj.t_query(end) - t_now) / (t_now - t_now_last_rec);
                        fprintf("%s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), Cr = (%s),\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), scheme_name, time_remained, join(string(size(X)), ", "), join(string(Cr), ", "), t_now, obj.t_query(t_query_next_idx), obj.t_query(end));
                        time_start = tic;
                        t_now_last_rec = t_now;
                    end

                    % add halo points, advance one time step
                    if scheme_name == "time_leap_frog_RA_space_central_diff"
                        if string(obj.bndry_func(0, 0)) == "periodic"
                            halo_left = f_last(end-1);
                            halo_right = f_last(2);
                        else
                            halo_left = obj.bndry_func(x_grid{1}(1) - obj.delta_x(1), t_last);
                            halo_right = obj.bndry_func(x_grid{1}(end) + obj.delta_x(1), t_last);
                        end
                        f_last_expanded = [halo_left; f_last; halo_right];
                        f_now = f_llast  ...
                                + Cr(1)*( + f_last_expanded(length(halo_left) : -1 + end-length(halo_right)) ...
                                          - f_last_expanded(2 + length(halo_left) : 1 + end-length(halo_right)) ...
                                        );
                        
                        % RA time filter
                        d_f = ra_opts.nu / 2 * (f_now - 2 * f_last + f_llast);
                        f_last = f_last + ra_opts.alpha * d_f;
                        f_now = f_now + (ra_opts.alpha - 1) * d_f;
                    end

                    % store solution at queried time.
                    if ~((obj.t_query(t_query_next_idx) - t_llast)*(obj.t_query(t_query_next_idx) - t_last) > 0)
                        if abs(obj.t_query(t_query_next_idx) - t_last) < abs(obj.t_query(t_query_next_idx) - t_llast)
                            t_list(t_query_next_idx) = t_last;
                            f_list{t_query_next_idx} = f_last;
                        else
                            t_list(t_query_next_idx) = t_llast;
                            f_list{t_query_next_idx} = f_llast;
                        end
                        t_query_next_idx = t_query_next_idx + 1;
                    end

                    % prepare for next time step
                    t_last = t_now;
                    f_llast = f_last; 
                    f_last = f_now;
                end % end of `while t_now < obj.t_query(end)`

                fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                return
            end % end of 

            warning("Invalid scheme_name. Please choose one of the following options: %s", join(AdvecProb2D.SCHEME_NAME, ", "));
        end % end of `function [f_list, x_grid, t_list] = solve_1D(obj, scheme_name)`

        function [f_list, x_grid, t_list] = solve_2D(obj, scheme_name)
            arguments (Input)
                obj
                scheme_name = "time_forward_Euler_space_upwind";
            end

            f_list = cell(size(obj.t_query));
            t_list = nan(size(obj.t_query));
            x_grid = cellfun(@(range_1D, delta_1D) (range_1D(1): delta_1D: range_1D(2)).', obj.x_range, num2cell(obj.delta_x), UniformOutput=false);
            if ismember(scheme_name, ["time_forward_Euler_space_upwind", "time_forward_Euler_space_CTU"])
                obj.scheme_name = scheme_name;

                Cr = abs(obj.velocity) * obj.delta_t ./ obj.delta_x;
                [X, Y] = ndgrid(x_grid{:});
                f_last = obj.init_func(X, Y);
                t_last = obj.t_start;
                t_now = t_last;
                t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

                % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
                if isempty(t_query_next_idx)
                    t_list(1) = t_last;
                    f_list{1} = f_last;

                    fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                    return
                end
                
                % 用于记录进度
                REC_TIME_INT = 5; % seconds
                time_start = tic;
                t_now_last_rec = t_last;

                while t_now < obj.t_query(end)
                    t_now = t_last + obj.delta_t;

                    % 输出进度
                    
                    if toc(time_start) > REC_TIME_INT
                        time_remained = REC_TIME_INT * (obj.t_query(end) - t_now) / (t_now - t_now_last_rec);
                        fprintf("%s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), Cr = (%s),\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), scheme_name, time_remained, join(string(size(X)), ", "), join(string(Cr), ", "), t_now, obj.t_query(t_query_next_idx), obj.t_query(end));
                        time_start = tic;
                        t_now_last_rec = t_now;
                    end
                    
                    % add halo points, advance one time step
                    if ismember(scheme_name, ["time_forward_Euler_space_upwind", "time_forward_Euler_space_CTU"])
                        if obj.bndry_func(0, 0, 0) == "periodic"
                            halo_left = [f_last(end-1, end-1); f_last(:, end-1); f_last(2, end-1)];
                            halo_right = [f_last(end-1, 2); f_last(:, 2); f_last(2, 2)];
                            halo_upper = f_last(end-1, :);
                            halo_lower = f_last(2, :);
                        else
                            halo_left = obj.bndry_func([x_grid{1}(1) - obj.delta_x(1); x_grid{1}; x_grid{1}(end) + obj.delta_x(1)], ones(2 + length(x_grid{1}), 1)*(x_grid{2}(1) - obj.delta_x(2)), t_last);
                            halo_right = obj.bndry_func([x_grid{1}(1) - obj.delta_x(1); x_grid{1}; x_grid{1}(end) + obj.delta_x(1)], ones(2 + length(x_grid{1}), 1)*(x_grid{2}(end) + obj.delta_x(2)), t_last);
                            halo_upper = obj.bndry_func(ones(1, length(x_grid{2}))*(x_grid{1}(1) - obj.delta_x(1)), x_grid{2}.', t_last);
                            halo_lower = obj.bndry_func(ones(1, length(x_grid{2}))*(x_grid{1}(end) + obj.delta_x(1)), x_grid{2}.', t_last);
                        end
                        f_last_expanded = [halo_left, [halo_upper; f_last; halo_lower], halo_right];
                        f_now = (1 - Cr(1) - Cr(2)) * f_last_expanded(2:end-1, 2:end-1) ...
                                + Cr(1)*(obj.velocity(1) > 0)*f_last_expanded(1:end-2, 2:end-1) ...
                                + Cr(1)*(obj.velocity(1) < 0)*f_last_expanded(3:end, 2:end-1) ...
                                + Cr(2)*(obj.velocity(2) > 0)*f_last_expanded(2:end-1, 1:end-2) ...
                                + Cr(2)*(obj.velocity(2) < 0)*f_last_expanded(2:end-1, 3:end);
                        if scheme_name == "time_forward_Euler_space_CTU"
                            f_now = f_now + prod(Cr) * f_last_expanded(2:end-1, 2:end-1) ...
                                    - prod(Cr) * (obj.velocity(1) > 0)*f_last_expanded(1:end-2, 2:end-1) ...
                                    - prod(Cr) * (obj.velocity(1) < 0)*f_last_expanded(3:end, 2:end-1) ...
                                    - prod(Cr) * (obj.velocity(2) > 0)*f_last_expanded(2:end-1, 1:end-2) ...
                                    - prod(Cr) * (obj.velocity(2) < 0)*f_last_expanded(2:end-1, 3:end) ...
                                    + prod(Cr) * (obj.velocity(1) > 0 & obj.velocity(2) > 0)*f_last_expanded(1:end-2, 1:end-2) ...
                                    + prod(Cr) * (obj.velocity(1) > 0 & obj.velocity(2) < 0)*f_last_expanded(1:end-2, 3:end) ...
                                    + prod(Cr) * (obj.velocity(1) < 0 & obj.velocity(2) > 0)*f_last_expanded(3:end, 1:end-2) ...
                                    + prod(Cr) * (obj.velocity(1) < 0 & obj.velocity(2) < 0)*f_last_expanded(3:end, 3:end);
                        end
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
                end % end of `while t_now < obj.t_query(end)`

                fprintf("%s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), scheme_name, join(string(t_list), ", "));
                return
            end

            warning("Invalid scheme_name. Please choose one of the following options: %s", join(AdvecProb2D.SCHEME_NAME, ", "));
        end % end of `function [f_list, x_grid, t_list] = solve_2D(obj, scheme_name)`
    end % end of `methods (Access=private)`
end % end of `classdef AdvecProb2D < handle`
