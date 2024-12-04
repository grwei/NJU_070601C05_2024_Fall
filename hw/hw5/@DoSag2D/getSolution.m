function [B_list, O_list, x_grid, y_grid, t_list, scheme_name] = getSolution(obj)
    if isempty(obj.f_list)
        error("Problem not solved yet. Please call obj.solve(scheme_name) before calling get_solution().")
    end

    B_list = arrayfun(@(t_ind) obj.f_list{t_ind}{1}, 1:length(obj.f_list), UniformOutput=false);
    O_list = arrayfun(@(t_ind) obj.f_list{t_ind}{2}, 1:length(obj.f_list), UniformOutput=false);
    x_grid = obj.x_grid;
    y_grid = obj.y_grid;
    t_list = obj.t_list;
    scheme_name = obj.scheme_name;
end
