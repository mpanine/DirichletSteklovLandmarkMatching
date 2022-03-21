clear all;
close all;

data_dir_glob = './output/007_lmks/';

data_dir = './data/FAUST/';
base_pair_file = [data_dir 'test_pairs.txt'];
pair_list = {};
fid = fopen(base_pair_file); % open the file
while ~feof(fid) % feof(fid) is true when the file ends
      textLineEntry = fgetl(fid); % read one line
      [pair,comma] = split(textLineEntry,'	');
      pair_list{end+1} = pair;
end
fclose(fid); % close the file

r_ls = [0.25,0.5,0.75];
%% Aggregate results for each # lmks and each method
thresh = 0:0.001:0.50; % evaluation thresholds

mean_curves = NaN*ones(length(r_ls),length(thresh));
mean_errs = NaN*ones(length(r_ls),1);
mean_times = NaN*ones(length(r_ls),1);

for r_idx=1:length(r_ls)
    r = r_ls(r_idx);
    
    all_curves = NaN*ones(length(pair_list),length(thresh));
    all_err = NaN*ones(length(pair_list),1);
    all_times = NaN*ones(length(pair_list),1);
    
    for pair_id=1:length(pair_list)
        p = pair_list{pair_id};
        s1 = p{1};
        s2 = p{2};
        %%
        data_dir = sprintf('%s/R_%.2f/',data_dir_glob,r);

        %%

        filename_str = [sprintf('%s/%s_%s_',data_dir,s1,s2) '%s.mat'];

        load(sprintf(filename_str,'curve'));
        load(sprintf(filename_str,'err'));
        load(sprintf(filename_str,'time'));

        all_curves(pair_id,:) = curve_ZO;
        all_err(pair_id) = mean(err_ZO);
        all_times(pair_id) = computation_time;
    end
    mean_curves(r_idx,:) = mean(all_curves,1);
    mean_errs(r_idx) = mean(all_err)
    mean_times(r_idx) = mean(all_times)
end
%%
f = figure;
plot(squeeze(r_ls),squeeze(mean_errs),'-r','lineWidth',2,...
                        'marker','o',...
                        'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r');
grid on;
xlim([min(r_ls),max(r_ls)]);
xticks(r_ls);
ylim([0.014,max(mean_errs)]);
xlabel('Mesh Reduction Factor');
ylabel('Mean Geodesic Error');
