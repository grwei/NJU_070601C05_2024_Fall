function params_struct = prepareParams(t_start, range_x, range_y, velocity, diffu_coeff, K_r, Q_B_loc, Q_B, K_a, O_s, init_func, bndry_func_type_x0, bndry_func_x0, bndry_func_type_x1, bndry_func_x1, bndry_func_type_y0, bndry_func_y0, bndry_func_type_y1, bndry_func_y1, delta_t, delta_x, delta_y, t_query)
    arguments
        t_start = 0;
        range_x (1, 2) = [0, 150];
        range_y (1, 2) = [0, 30];
        velocity (1, 2) = [.4, 0];
        diffu_coeff (1, 2) = [.5, .5];
        K_r {isscalar} = .01;
        Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
        Q_B = [70; 0];    % BOD 点源强度
        K_a {isscalar} = .02;
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
        delta_t {isscalar} = .01;
        delta_x {isscalar} = .2;
        delta_y {isscalar} = .2;
        t_query {isvector} = [0, 1, 2, +Inf];
    end

    % 求解区域
    params_struct.t_start = t_start;
    params_struct.range_x = sort(range_x, "ascend");
    params_struct.range_y = sort(range_y, "ascend"); 

    params_struct.velocity = velocity;
    params_struct.diffu_coeff = diffu_coeff;
    params_struct.K_r = K_r;
    params_struct.Q_B_loc = Q_B_loc;
    params_struct.Q_B = Q_B;
    params_struct.K_a = K_a;
    params_struct.O_s = O_s;

    params_struct.init_func = init_func;
    params_struct.bndry_func_type_x0 = bndry_func_type_x0;
    params_struct.bndry_func_x0 = bndry_func_x0;
    params_struct.bndry_func_type_x1 = bndry_func_type_x1;
    params_struct.bndry_func_x1 = bndry_func_x1;
    params_struct.bndry_func_type_y0 = bndry_func_type_y0;
    params_struct.bndry_func_y0 = bndry_func_y0;
    params_struct.bndry_func_type_y1 = bndry_func_type_y1;
    params_struct.bndry_func_y1 = bndry_func_y1;

    params_struct.delta_t = delta_t;
    params_struct.delta_x = delta_x;
    params_struct.delta_y = delta_y;
    params_struct.t_query = unique(t_query(t_query >= t_start), "sorted");
end
