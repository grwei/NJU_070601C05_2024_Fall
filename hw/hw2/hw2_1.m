%% hw2_1.m
% Description: 用 Thomas algorithm (追赶法) 求解系数矩阵为三对角的线性方程组
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 17, 2024
% Last modified: Oct. 18, 2024
%

%%

clear; clc; close all

%% Problem 1(a)

%%% Constant

N = 20;
phi_0 = 2;
phi_N = 3;

%%% solve

tri_diag_coeff_mat = ones(N-1, 1) * [1, -2, 1];
const_vec = zeros(N-1, 1);
const_vec([1, end]) = -[phi_0, phi_N];

solver = TriDiagSolver(tri_diag_coeff_mat, const_vec);
phi_numeric = [phi_0; solver.solve(); phi_N];
solver.get_cpu_time()

%%% figure

t_fig = figure(Name="fig_1_1_laplace_eqn_thomas_results");

% set figure size
UNIT_ORIGINAL = t_fig.Units;
t_fig.Units = "centimeters";
t_fig.Position = [3, 3, 12, 12];
t_fig.Units = UNIT_ORIGINAL;

% create figure
t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="compact");
t_axes = nexttile(t_TCL, 1);

hold(t_axes, "on")
plot(t_axes, linspace(0, 1, N+1).', phi_numeric, LineWidth=1.5, DisplayName="(numeric) $\phi_{xx} = 0, \, \phi(0) = 2, \, \phi(1) = 3$")
plot(t_axes, [0, 1], [2, 3], '--k', LineWidth=1.5, DisplayName="(exact) $\phi = \phi(0) + (\phi(1) - \phi(0)) x$")
set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", XLimitMethod="tight")

legend(t_axes, Location="best", Interpreter="latex", FontSize=10.5, Box="off")
xlabel(t_axes, "$x$", Interpreter="latex", FontSize=10.5);
ylabel(t_axes, "$\phi$", Interpreter="latex", FontSize=10.5);

print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")

%% Problem 1(b)

% for N_2_exp_max = [12, 16]
%     hw2_1b(N_2_exp_max, true)
% end

for N_2_exp_max = [11, 20]
    hw2_1b(N_2_exp_max, true)
end

%%

function hw2_1b(N_2_exp_max, flag_save_fig)
    arguments
        N_2_exp_max = 11;
        flag_save_fig = false;
    end
    N = floor(2.^(5:0.25:N_2_exp_max));
    t_thomas = nan(size(N));
    t_sparse = nan(size(N));
    t_full = nan(size(N));
    t_randi = nan(size(N));
    for i = 1:length(N)
        [t_thomas(i), t_sparse(i), t_full(i), t_randi(i)] = hw2_1_2(N(i));
    end
    
    %%% linear fit
    
    lin_fit_thomas = fitlm(log2(N), log2(t_thomas)).Coefficients;
    lin_fit_sparse = fitlm(log2(N), log2(t_sparse)).Coefficients;
    lin_fit_full = fitlm(log2(N), log2(t_full)).Coefficients;
    lin_fit_randi = fitlm(log2(N), log2(t_randi)).Coefficients;
    
    %%% figure
    
    t_fig = figure(Name="fig_1_2_linear_eqns_solvers_comparison_" + string(N_2_exp_max));
    
    % set figure size
    UNIT_ORIGINAL = t_fig.Units;
    t_fig.Units = "centimeters";
    t_fig.Position = [3, 3, 12, 12];
    t_fig.Units = UNIT_ORIGINAL;
    
    % create figure
    t_TCL = tiledlayout(t_fig, 1, 1, TileSpacing="compact", Padding="compact");
    t_axes = nexttile(t_TCL, 1);
    
    hold(t_axes, "on")
    plot(t_axes, log2(N), log2(t_thomas), LineWidth=1.5, DisplayName=sprintf("thomas\\_algo $(%.2f \\pm %.2f)$", lin_fit_thomas.Estimate(2), lin_fit_thomas.SE(2)))
    plot(t_axes, log2(N), log2(t_sparse), LineWidth=1.5, DisplayName=sprintf("mldivide\\_sparse $(%.2f \\pm %.2f)$", lin_fit_sparse.Estimate(2), lin_fit_sparse.SE(2)))
    plot(t_axes, log2(N), log2(t_full), '--', LineWidth=1.5, DisplayName=sprintf("mldivide\\_full $(%.2f \\pm %.2f)$", lin_fit_full.Estimate(2), lin_fit_full.SE(2)))
    plot(t_axes, log2(N), log2(t_randi), LineWidth=1.5, DisplayName=sprintf("mldivide\\_randi $(%.2f \\pm %.2f)$", lin_fit_randi.Estimate(2), lin_fit_randi.SE(2)))
    set(t_axes, FontName="Times New Roman", FontSize=10.5, Box="on", TickLabelInterpreter="latex", XLimitMethod="tight")
    
    legend(t_axes, Location="best", Interpreter="latex", FontSize=10.5, Box="off")
    xlabel(t_axes, "$\log_2{N}$", Interpreter="latex", FontSize=10.5);
    ylabel(t_axes, "$\log_2{T}$", Interpreter="latex", FontSize=10.5);
    
    print(t_fig, ".\fig\" + t_fig.Name + ".svg", "-dsvg")                                                                                                                      
end

%%

function [t_thomas, t_sparse, t_full, t_randi] = hw2_1_2(N)
    arguments
        N = 2^10
    end

    vec_rand = randi(100, N, 4);
    t_thomas = timeit(@() thomas());
    mat_coeff_sparse = spdiags(vec_rand(:, 1:3), -1:1, N, N);
    t_sparse = timeit(@() mat_coeff_sparse \ vec_rand(:, 4));

    if N <= 2^15
        mat_coeff_full = full(mat_coeff_sparse);
        wrap_full = @() mat_coeff_full \ vec_rand(:, 4);
        t_full = timeit(wrap_full);
    else
        t_full = NaN;
    end
    
    if N <= 2^13
        mat_coeff_randi = randi(100, N);
        wrap_randi = @() mat_coeff_randi \ vec_rand(:, 4);
        t_randi = timeit(wrap_randi);
    else
        t_randi = NaN;
    end

    function [solution] = thomas()
        solver = TriDiagSolver([[NaN; vec_rand(1:end-1, 1)], vec_rand(:, 2), [vec_rand(2:end, 3); NaN]], vec_rand(:, 4));
        solution = solver.solve();
    end
end
