function [f_list, x_grid, y_grid, t_list] = solve2D(obj, scheme_name)
    arguments (Input)
        obj
        scheme_name = "time_forward_Euler_space_central";
    end

    x_grid = obj.x_grid;
    y_grid = obj.y_grid;
    [X, Y] = ndgrid(obj.x_grid, obj.y_grid);
    f_list = cell(size(obj.t_query));
    t_list = nan(size(obj.t_query));
    
    if ismember(scheme_name, "time_forward_Euler_space_central")
        obj.scheme_name = scheme_name;

        Cr = obj.velocity * obj.delta_t ./ [obj.delta_x, obj.delta_y];
        r = obj.diffu_coeff * obj.delta_t ./ [obj.delta_x, obj.delta_y].^2;
        
        f_last = cellfun(@(i_func) i_func(X, Y), obj.init_func, UniformOutput=false);
        % f_last{1} = obj.init_func{1}(X, Y); % BOD
        % f_last{2} = obj.init_func{2}(X, Y); % DO
        f_now = f_last;
        t_last = obj.t_start;
        t_query_next_idx = find(obj.t_query > obj.t_start, 1, "first");

        % 若请求时刻不晚于初始时刻, 则无需计算, 直接返回初值
        if isempty(t_query_next_idx) || t_query_next_idx > 1
            t_list(1) = t_last;
            f_list{1} = f_last;

            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tInitial value added.\n\tt_start = %g\n", datetime('now'), obj.desc_str, scheme_name, t_last);
            fclose(LOG_FILE_ID);
            fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tInitial value added.\n\tt_start = %g\n", datetime('now'), obj.desc_str, scheme_name, t_last);
        end

        %%% 若请求 +Inf 时刻, 则创建稳态解, obj.f_steady
        if any(obj.t_query == +Inf)
            f_steady = obj.solve2DSteady(scheme_name);
            t_list(obj.t_query == +Inf) = +Inf;
            f_list{obj.t_query == +Inf} = f_steady;

            LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
            fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tSteady state solved.\n\tt_steady = +Inf, DO_min (sat) = %.3g (%.3g),\n", datetime('now'), obj.desc_str, scheme_name, min(f_steady{2}, [], "all"), obj.O_s);
            fclose(LOG_FILE_ID);
            fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tSteady state solved.\n\tt_steady = +Inf, DO_min (sat) = %.3g (%.3g),\n", datetime('now'), obj.desc_str, scheme_name, min(f_steady{2}, [], "all"), obj.O_s);
        end

        %%% 下面求解有限时刻的解.

        T_END = obj.t_query(obj.t_query < Inf);
        T_END = T_END(end);
        
        % 用于记录进度
        REC_TIME_INT = 5; % seconds
        time_start = tic;
        t_now_last_rec = t_last;

        % 测试
        if ~isfolder(".\log\")
            mkdir(".\log\")
        end
        thres = [.5, (obj.O_s - 6) / (obj.O_s - min(f_steady{2}, [], "all")), .99, .999];
        flag_test = true(size(thres));

        while t_last < T_END
            t_now = t_last + obj.delta_t;

            % 输出进度
            if toc(time_start) > REC_TIME_INT
                test_prog = (obj.O_s - min(f_last{2}, [], "all")) / (obj.O_s - min(f_steady{2}, [], "all"));
                time_remained = REC_TIME_INT * (T_END - t_now) / (t_now - t_now_last_rec);
                fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\ttime_remained = %.2e,\n\tspace_dim = (%s), Cr = (%s), r = (%s), test_prog = %.2f %%,\n\tt_now = %.2e, t_query_next = %.2e, t_query_end = %.2e\n", datetime('now'), obj.desc_str, scheme_name, time_remained, join(string(size(X)), ", "), join(string(Cr), ", "), join(string(r), ", "), test_prog*100, t_now, obj.t_query(t_query_next_idx), T_END);
                time_start = tic;
                t_now_last_rec = t_now;
            end
            
            % add halo points, advance one time step
            if scheme_name == "time_forward_Euler_space_central"
                switch obj.bndry_func_type_x0
                    case "Dirichlet"
                        halo_upper = cellfun(@(bndry_func, f) 2*bndry_func(obj.y_grid.', t_last) - f(2, :), obj.bndry_func_x0, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_upper = cellfun(@(bndry_func, f) f(2, :) - 2*obj.delta_x*bndry_func(obj.y_grid.', t_last), obj.bndry_func_x0, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_upper = cellfun(@(f) f(end-1, :), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_x1
                    case "Dirichlet"
                        halo_lower = cellfun(@(bndry_func, f) 2*bndry_func(obj.y_grid.', t_last) - f(end-1, :), obj.bndry_func_x1, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_lower = cellfun(@(bndry_func, f) f(end-1, :) + 2*obj.delta_x*bndry_func(obj.y_grid.', t_last), obj.bndry_func_x1, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_lower = cellfun(@(f) f(2, :), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_y0
                    case "Dirichlet"
                        halo_left = cellfun(@(bndry_func, f) 2*bndry_func(obj.x_grid, t_last) - f(:, 2), obj.bndry_func_y0, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_left = cellfun(@(bndry_func, f) f(:, 2) - 2*obj.delta_y*bndry_func(obj.x_grid, t_last), obj.bndry_func_y0, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_left = cellfun(@(f) f(:, end-1), f_last, UniformOutput=false);
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end

                switch obj.bndry_func_type_y1
                    case "Dirichlet"
                        halo_right = cellfun(@(bndry_func, f) 2*bndry_func(obj.x_grid, t_last) - f(:, end-1), obj.bndry_func_y1, f_last, UniformOutput=false);
                    case "Neumann"
                        halo_right = cellfun(@(bndry_func, f) f(:, end-1) + 2*obj.delta_y*bndry_func(obj.x_grid, t_last), obj.bndry_func_y1, f_last, UniformOutput=false);
                    case "Periodic"
                        halo_right = cellfun(@(f) f(:, 2), f_last, UniformOutput=false) ;
                    otherwise
                        error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
                end
                
                f_last_expanded = cellfun(@(f_l, h_left, h_right, h_upper, h_lower) [[NaN; h_left; NaN], [h_upper; f_l; h_lower], [NaN; h_right; NaN]], f_last, halo_left, halo_right, halo_upper, halo_lower, UniformOutput=false);
                f_now{1} = (1 - 2*r(1) - 2*r(2) - obj.delta_t*obj.K_r) * f_last{1} ...
                         + (r(1) - Cr(1)/2) * f_last_expanded{1}(3:end, 2:end-1) ...
                         + (r(1) + Cr(1)/2) * f_last_expanded{1}(1:end-2, 2:end-1) ...
                         + (r(2) - Cr(2)/2) * f_last_expanded{1}(2:end-1, 3:end) ...
                         + (r(2) + Cr(2)/2) * f_last_expanded{1}(2:end-1, 1:end-2) ...
                         + obj.delta_t * obj.B_src;
                f_now{2} = (1 - 2*r(1) - 2*r(2) - obj.delta_t*obj.K_a)*f_last{2} ...
                + (r(1) - Cr(1)/2) * f_last_expanded{2}(3:end, 2:end-1) ...
                + (r(1) + Cr(1)/2) * f_last_expanded{2}(1:end-2, 2:end-1) ...
                + (r(2) - Cr(2)/2) * f_last_expanded{2}(2:end-1, 3:end) ...
                + (r(2) + Cr(2)/2) * f_last_expanded{2}(2:end-1, 1:end-2) ...
                + obj.delta_t * (obj.K_a * obj.O_s - obj.K_r * f_last{1});
            else
                error("格式 %s 正在开发", scheme_name)
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

            % 测试
            test_prog = (obj.O_s - min(f_last{2}, [], "all")) / (obj.O_s - min(f_steady{2}, [], "all"));
            for test_ind = 1:length(thres)
                if flag_test(test_ind) && test_prog > thres(test_ind)
                    flag_test(test_ind) = false;
                    LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
                    fprintf(LOG_FILE_ID, "%s\n\tt = %g, m = %.3g (m_s = %.3g), prog = %.3f %%\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), t_last, min(f_last{2}, [], "all"), min(f_steady{2}, [], "all"), 100*test_prog);
                    fprintf("%s\n\tt = %g, m = %.3g (m_s = %.3g), prog = %.3f %%\n", string(datetime("now"), "yyyy-MM-dd HH:mm:ss", "en_US"), t_last, min(f_last{2}, [], "all"), min(f_steady{2}, [], "all"), 100*test_prog);
                    fclose(LOG_FILE_ID);
                end
            end
        end % end of `while t_now < obj.t_query(end)`

        LOG_FILE_ID = fopen(obj.LOG_FILE_PATH,'a');
        fprintf(LOG_FILE_ID, "%s\n\tdescription: %s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), obj.desc_str, scheme_name, join(string(t_list), ", "));
        fclose(LOG_FILE_ID);
        fprintf("%s\n\tdescription: %s\n\tscheme_name = %s,\n\tFinished.\n\tt_list = [%s]\n", datetime('now'), obj.desc_str, scheme_name, join(string(t_list), ", "));
        return
    end % end of `if ismember(scheme_name, "time_forward_Euler_space_central")`

    error("Invalid scheme_name. Please choose one of the following options: %s", join(DoSag2D.SCHEME_NAME, ", "));
end % end of `function [f_list, x_grid, t_list] = solve_2D(obj, scheme_name)`
