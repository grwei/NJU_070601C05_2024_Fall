function f_steady = solve2DSteady(obj, scheme_name)
    arguments (Input)
        obj
        scheme_name = "time_forward_Euler_space_central";
    end

    if ~isempty(obj.f_steady)
        f_steady = obj.f_steady;
        return
    end
    f_steady = cell([2, 1]);

    if scheme_name == "time_forward_Euler_space_central"
        Cr = obj.velocity * obj.delta_t ./ [obj.delta_x, obj.delta_y];
        r = obj.diffu_coeff * obj.delta_t ./ [obj.delta_x, obj.delta_y].^2;

        %%% 1. 先求 B_steady
        
        % 1.1 创建方程组的系数矩阵 (sparse) 和常向量 (full)

        sz = [length(obj.x_grid), length(obj.y_grid)];
        c_vec{2} = nan(prod(sz), 1);
        c_vec{1} = c_vec{2};
        i_list{2} = nan(prod(sz - 2)*5 + sum(sz - 2)*6 + 8, 1);
        i_list{1} = i_list{2};
        j_list = i_list;
        v_list = i_list;
        cnt = 0;

        % 1.1.1 建立 halo_points
        switch obj.bndry_func_type_x0
            case "Dirichlet"
                halo_upper = {2, [2, -1]}; % halo point 是 边界函数 以及 法向格点值 的 线性组合, 记作 {边界函数 系数, [格点 1 索引, 系数 1; 格点 2 索引, 系数 2]}
            case "Neumann"
                halo_upper = {-2*obj.delta_x, [2, 1]};
            case "Periodic"
                % halo_upper = {0, [length(obj.x_grid)-1, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_x1
            case "Dirichlet"
                halo_lower = {2, [length(obj.x_grid)-1, -1]};
            case "Neumann"
                halo_lower = {2*obj.delta_x, [length(obj.x_grid)-1, 1]};
            case "Periodic"
                % halo_lower = {0, [2, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_y0
            case "Dirichlet"
                halo_left = {2, [2, -1]};
            case "Neumann"
                halo_left = {-2*obj.delta_y, [2, 1]};
            case "Periodic"
                % halo_left = {0, [length(obj.y_grid)-1, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        switch obj.bndry_func_type_y1
            case "Dirichlet"
                halo_right = {2, [length(obj.y_grid)-1, -1]};
            case "Neumann"
                halo_right = {2*obj.delta_y, [length(obj.y_grid)-1, 1]};
            case "Periodic"
                % halo_right = {0, [2, 1]};
                error("周期边界条件下的稳态解程序未完成")
            otherwise
                error("不支持的边界条件类型: %s. 请选择一种: [%s]", obj.bndry_func_type_x0, join(["Dirichlet", "Neumann", "Periodic"], ", "))
        end

        % 1.1.2 建立 BOD 方程 (以后, 计划把 halo points 也写进方程的未知数, 便于统一各类边界条件的实现)
        for ind_x = 1:sz(1)
            for ind_y = 1:sz(2)
                ind_x_guess = ind_x + [0; 1; -1; 0; 0]; % 以后: 应加入用于表示四个 halo 的格点
                ind_y_guess = ind_y + [0; 0; 0; 1; -1];
                flag_invalid = ind_x_guess < 1 | ind_x_guess > sz(1) | ind_y_guess < 1 | ind_y_guess > sz(2);
                
                ind_x_valid = ind_x_guess(~flag_invalid);
                ind_y_valid = ind_y_guess(~flag_invalid);
                ind_lin = sub2ind(sz, ind_x_valid, ind_y_valid);
                c_vec{1}(ind_lin(1)) = obj.delta_t * obj.B_src(ind_x, ind_y) ...
                    + (ind_x == 1)*(r(1) + Cr(1)/2)*halo_upper{1}*obj.bndry_func_x0{1}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_x == sz(1))*(r(1) - Cr(1)/2)*halo_lower{1}*obj.bndry_func_x1{1}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_y == 1)*(r(2) + Cr(2)/2)*halo_left{1}*obj.bndry_func_y0{1}(obj.x_grid(ind_x), +Inf) ...
                    + (ind_y == sz(2))*(r(2) - Cr(2)/2)*halo_right{1}*obj.bndry_func_y1{1}(obj.x_grid(ind_x), +Inf);
                
                nnz_ = length(ind_lin);
                i_list{1}(cnt + 1 : cnt + nnz_) = ind_lin(1);
                j_list{1}(cnt + 1 : cnt + nnz_) = ind_lin;
                v_list_ = [ind_x_valid == ind_x & ind_y_valid == ind_y, ind_x_valid == ind_x+1, ind_x_valid == ind_x-1, ind_y_valid == ind_y+1, ind_y_valid == ind_y-1] * [ ...
                    2*sum(r) + obj.delta_t * obj.K_r; % B_{i,j}
                    -r(1) + Cr(1)/2; % B_{i+1, j}
                    -r(1) - Cr(1)/2; % B_{i-1, j}
                    -r(2) + Cr(2)/2; % B_{i, j+1}
                    -r(2) - Cr(2)/2  % B_{i, j-1}
                ];
                if ind_x == 1
                    % 以后: 用循环, 考虑用于表示 halo 的所有格点
                    v_list_(ind_x_valid == halo_upper{2}(1,1)) = v_list_(ind_x_valid == ind_x + 1) - (r(1) + Cr(1)/2)*halo_upper{2}(1,2);
                end
                if ind_x == sz(1)
                    v_list_(ind_x_valid == halo_lower{2}(1,1)) = v_list_(ind_x_valid == ind_x - 1) - (r(1) - Cr(1)/2)*halo_lower{2}(1,2);
                end
                if ind_y == 1
                    v_list_(ind_y_valid == halo_left{2}(1,1)) = v_list_(ind_y_valid == ind_y + 1) - (r(2) + Cr(2)/2)*halo_left{2}(1,2);
                end
                if ind_y == sz(2)
                    v_list_(ind_y_valid == halo_right{2}(1,1)) = v_list_(ind_y_valid == ind_y - 1) - (r(2) - Cr(2)/2)*halo_right{2}(1,2);
                end
                v_list{1}(cnt + 1 : cnt + nnz_) = v_list_;

                cnt = cnt + nnz_;
            end
        end

        if cnt ~= length(v_list{1})
            error("这是不可能的, 除非程序有误")
        end

        % 1.1.3 求解 BOD 稳态方程

        B_steady = sparse(i_list{1}, j_list{1}, v_list{1}, prod(sz), prod(sz)) \ c_vec{1};
        f_steady{1} = reshape(B_steady, sz);

        %%% 2. 再求 O_steady

        i_list{2} = i_list{1};
        j_list{2} = j_list{1};

        cnt = 0;
        for ind_x = 1:sz(1)
            for ind_y = 1:sz(2)
                ind_x_guess = ind_x + [0; 1; -1; 0; 0]; % 以后: 应加入用于表示四个 halo 的格点
                ind_y_guess = ind_y + [0; 0; 0; 1; -1];
                flag_invalid = ind_x_guess < 1 | ind_x_guess > sz(1) | ind_y_guess < 1 | ind_y_guess > sz(2);
                
                ind_x_valid = ind_x_guess(~flag_invalid);
                ind_y_valid = ind_y_guess(~flag_invalid);
                ind_lin = sub2ind(sz, ind_x_valid, ind_y_valid);
                c_vec{2}(ind_lin(1)) = obj.delta_t * (obj.K_a * obj.O_s - obj.K_r * f_steady{1}(ind_x, ind_y)) ...
                    + (ind_x == 1)*(r(1) + Cr(1)/2)*halo_upper{1}*obj.bndry_func_x0{2}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_x == sz(1))*(r(1) - Cr(1)/2)*halo_lower{1}*obj.bndry_func_x1{2}(obj.y_grid(ind_y), +Inf) ...
                    + (ind_y == 1)*(r(2) + Cr(2)/2)*halo_left{1}*obj.bndry_func_y0{2}(obj.x_grid(ind_x), +Inf) ...
                    + (ind_y == sz(2))*(r(2) - Cr(2)/2)*halo_right{1}*obj.bndry_func_y1{2}(obj.x_grid(ind_x), +Inf);
                
                nnz_ = length(ind_lin);
                v_list_ = [ind_x_valid == ind_x & ind_y_valid == ind_y, ind_x_valid == ind_x+1, ind_x_valid == ind_x-1, ind_y_valid == ind_y+1, ind_y_valid == ind_y-1] * [ ...
                    2*sum(r) + obj.delta_t * obj.K_a; % O_{i,j}
                    -r(1) + Cr(1)/2; % O_{i+1, j}
                    -r(1) - Cr(1)/2; % O_{i-1, j}
                    -r(2) + Cr(2)/2; % O_{i, j+1}
                    -r(2) - Cr(2)/2  % O_{i, j-1}
                ];
                if ind_x == 1
                    % 以后: 用循环, 考虑用于表示 halo 的所有格点
                    v_list_(ind_x_valid == halo_upper{2}(1,1)) = v_list_(ind_x_valid == ind_x + 1) - (r(1) + Cr(1)/2)*halo_upper{2}(1,2);
                end
                if ind_x == sz(1)
                    v_list_(ind_x_valid == halo_lower{2}(1,1)) = v_list_(ind_x_valid == ind_x - 1) - (r(1) - Cr(1)/2)*halo_lower{2}(1,2);
                end
                if ind_y == 1
                    v_list_(ind_y_valid == halo_left{2}(1,1)) = v_list_(ind_y_valid == ind_y + 1) - (r(2) + Cr(2)/2)*halo_left{2}(1,2);
                end
                if ind_y == sz(2)
                    v_list_(ind_y_valid == halo_right{2}(1,1)) = v_list_(ind_y_valid == ind_y - 1) - (r(2) - Cr(2)/2)*halo_right{2}(1,2);
                end
                v_list{2}(cnt + 1 : cnt + nnz_) = v_list_;

                cnt = cnt + nnz_;
            end
        end

        if cnt ~= length(v_list{2})
            error("这是不可能的, 除非程序有误")
        end

        O_steady = sparse(i_list{2}, j_list{2}, v_list{2}, prod(sz), prod(sz)) \ c_vec{2};
        f_steady{2} = reshape(O_steady, sz);

        %%% 3.

        obj.f_steady = f_steady;
    else
        error("scheme_name: %s 的稳态解算法未实现", scheme_name)
    end
end
