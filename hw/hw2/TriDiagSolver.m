%% TriDiagSolver.m
% Description: 用 Thomas algorithm (追赶法) 求解系数矩阵为三对角的线性方程组
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 17, 2024
% Last modified: Oct. 18, 2024
%

%% Class definition

classdef TriDiagSolver < handle
    properties (Access=private)
        tri_diag_coeff_mat__ = [];
        const_vec__ = [];
        solut_vec__ = [];
        cpu_time__ = [];
    end

    methods (Access=public)
        function obj = TriDiagSolver(tri_diag_coeff_mat, const_vec)
            arguments
                tri_diag_coeff_mat = [NaN, 2, 1;
                                      1, 2, 1;
                                      1, 2, NaN];
                const_vec = [4; 8; 8];
            end
            obj.tri_diag_coeff_mat__ = tri_diag_coeff_mat;
            obj.const_vec__ = const_vec;
        end

        function tri_diag_coeff_mat = get_coeff_mat(obj)
            tri_diag_coeff_mat = obj.tri_diag_coeff_mat__; 
        end

        function const_vec = get_const_vec(obj)
            const_vec = obj.const_vec__; 
        end

        function cpu_time = get_cpu_time(obj)
            cpu_time = obj.cpu_time__;
        end

        function obj = reset(obj, tri_diag_coeff_mat, const_vec)
            arguments
                obj
                tri_diag_coeff_mat = [NaN, 2, 1;
                                      1, 2, 1;
                                      1, 2, NaN];
                const_vec = [4; 8; 8];
            end
            obj.clear();
            obj.tri_diag_coeff_mat__ = tri_diag_coeff_mat;
            obj.const_vec__ = const_vec;
        end

        function obj = clear(obj)
            obj.tri_diag_coeff_mat__ = [];
            obj.const_vec__ = [];
            obj.solut_vec__ = [];
            obj.cpu_time__ = [];
        end
        
        function solut_vec = solve(obj)
            if ~isempty(obj.solut_vec__)
                solut_vec = obj.solut_vec__;
                return
            end

            aug_mat = [obj.tri_diag_coeff_mat__, obj.const_vec__];
            tStart = cputime();
            for i = 2:size(aug_mat, 1)
                aug_mat(i, [1, 2, 4]) = aug_mat(i, [1, 2, 4]) + aug_mat(i-1, [2, 3, 4]) * (-aug_mat(i, 1) / aug_mat(i-1, 2)); 
            end
            
            obj.solut_vec__(length(obj.const_vec__), 1) = aug_mat(length(obj.const_vec__), 4) / aug_mat(length(obj.const_vec__), 2);
            for i = length(obj.const_vec__)-1 : -1 : 1
                obj.solut_vec__(i) = (aug_mat(i, 4) - aug_mat(i, 3) * obj.solut_vec__(i + 1)) / aug_mat(i, 2);
            end

            obj.cpu_time__ = cputime() - tStart;
            solut_vec = obj.solut_vec__;
        end
    end % end of method
end % end of classdef
