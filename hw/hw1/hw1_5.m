%% hw1_5.m
% Description: create figure of the results of modified wave number analysis
% Author: Guorui Wei (危国锐) (313017602@qq.com)
% Created at: Oct. 16, 2024
% Last modified:
%

mod_wave_func_handle = {
    @(x) 4 * sin(sqrt(x) / 2).^2;
    @(x) 6 - 18 ./ (3 + 2 * tan(sqrt(x) / 2).*2);
};
scheme_disp_name_ = {
    "(5.16) $D_2: \; y = 4 \sin^2 (\sqrt{x} / 2)$";
    "(5.22) $D_4: \; y = 6 - 18 / (3 + 2 \tan^2 (\sqrt{x} / 2))$";
};

hw1_5_obj = ModWaveNumAnalysis();
hw1_5_obj.set_mod_wave_func_handle(mod_wave_func_handle);
hw1_5_obj.set_scheme_disp_name(scheme_disp_name_);
hw1_5_obj.create_fig(true);
