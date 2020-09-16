
function sim1_wrapMMAD(mix_name_prefix, ref_name_prefix, dat_path, output_path)
addpath(genpath('/path/to/MMAD/'));
nMarker = 3;
nDataset = 3;
nGrid = 10;

for m = 1:nMarker
    for d = 1:nDataset
        for p = 1:nGrid
            mix_file = strcat(dat_path,'/', mix_name_prefix,'_D',num2str(d),'P',num2str(p),'.txt')
            ref_file = strcat(dat_path,'/',ref_name_prefix,'_D',num2str(m),'.txt')
            out_file = strcat(output_path,'/','Results_',mix_name_prefix,'_D',num2str(d),'P',num2str(p),'M',num2str(m),'.txt')
            HJver_MMAD_MarkerDeconv(mix_file, ref_file, out_file)
        end
    end
end

end
