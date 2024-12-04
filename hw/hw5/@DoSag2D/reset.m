function obj = reset(obj, params_struct)
    field_missing_flag = ~isfield(params_struct, DoSag2D.PARAMS_REQUIRED);
    if any(field_missing_flag)
        warning("Invalid params_struct. Params not changed. Missing field: %s", join(DoSag2D.PARAMS_REQUIRED(field_missing_flag), ", "));
        return
    end

    params_struct = DoSag2D.prepareParams(params_struct.t_start, params_struct.range_x, params_struct.range_y, params_struct.velocity, params_struct.diffu_coeff, params_struct.K_r, params_struct.Q_B_loc, params_struct.Q_B, params_struct.K_a, params_struct.O_s, params_struct.init_func, params_struct.bndry_func_type_x0, params_struct.bndry_func_x0, params_struct.bndry_func_type_x1, params_struct.bndry_func_x1, params_struct.bndry_func_type_y0, params_struct.bndry_func_y0, params_struct.bndry_func_type_y1, params_struct.bndry_func_y1, params_struct.delta_t, params_struct.delta_x, params_struct.delta_y, params_struct.t_query);
    for field_name = DoSag2D.PARAMS_REQUIRED
        obj.(field_name) = params_struct.(field_name);
    end

    obj.x_grid = (obj.range_x(1): obj.delta_x: obj.range_x(2)).';
    obj.y_grid = (obj.range_y(1): obj.delta_y: obj.range_y(2)).';
    obj.B_src = obj.bSrc();

    obj.t_list = [];
    obj.f_steady = [];
    obj.f_list = [];
end
