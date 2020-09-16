
function sim2_wrapMMAD(mix_name_prefix, ref_name_prefix, dat_path, output_path)
addpath(genpath('/path/to/MMAD/'));
nComp = 6;

for c = 1:nComp
    mix_file = strcat(dat_path,'/', mix_name_prefix,'_C',num2str(c),'.txt')
    ref_file = strcat(dat_path,'/',ref_name_prefix,'_C',num2str(c),'.txt')
    out_file = strcat(output_path,'/','Results_',mix_name_prefix,'_C',num2str(c),'.txt')
    HJver_MMAD_MarkerDeconv(mix_file, ref_file, out_file)

end

end
