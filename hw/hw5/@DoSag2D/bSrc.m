function B_src = bSrc(obj)
    % Q_B_loc = [10,15; 20,30]; % [x1,y1; x2, y2;] BOD 点源位置
    % Q_B = [70; 0];    % BOD 点源强度

    x_bin_edges = [obj.x_grid - obj.delta_x / 2; obj.x_grid(end) + obj.delta_x / 2];
    x_bin_ind = discretize(obj.Q_B_loc(:, 1), x_bin_edges);

    y_bin_edges = [obj.y_grid - obj.delta_y / 2; obj.y_grid(end) + obj.delta_y / 2];
    y_bin_ind = discretize(obj.Q_B_loc(:, 2), y_bin_edges);

    B_src = zeros([length(obj.x_grid), length(obj.y_grid)]);
    for src_ind = 1:length(obj.Q_B)
        B_src(x_bin_ind(src_ind), y_bin_ind(src_ind)) = B_src(x_bin_ind(src_ind), y_bin_ind(src_ind)) + obj.Q_B(src_ind) / obj.delta_x / obj.delta_y;
    end
end
