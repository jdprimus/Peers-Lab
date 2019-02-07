%% Create Nonlinear Data table
% Jeremy Primus
% February 6, 2019

%% description
% script to reproducibly create table of nonlinear RNA data from files
% provided by Wen Zhou


%%
% import and initialize
clear all
all_rpkm = readtable("C:\Users\Peers Lab\Desktop\Jeremy - Spring 2019\Synechocystis\data_rpkm_filter.csv");
nonlinear_list = readtable("C:\Users\Peers Lab\Desktop\Jeremy - Spring 2019\Synechocystis\nonlinear_list.csv");

% cross-reference lists
[~, ia, ~] = intersect(all_rpkm.Name, nonlinear_list.Var1);

% extract nonlinear data only
nonlinear_rpkm_allreps = all_rpkm(ia,:);

[~, ncols] = size(nonlinear_rpkm_allreps);

for i = 2:ncols
    short_headings{i,:} = nonlinear_rpkm_allreps.Properties.VariableNames{i}(1:6);
end

utimepoints = unique(short_headings(2:end),'stable');
for j = 1:length(utimepoints)
    nonlinear_rpkm_avgs(:,j) = mean(nonlinear_rpkm_allreps{:,strcmp(short_headings, utimepoints{j})},2);
end
nonlinear_rpkm_avgs = array2table(nonlinear_rpkm_avgs);
nonlinear_rpkm_avgs = [nonlinear_rpkm_allreps.Name nonlinear_rpkm_avgs];
nonlinear_rpkm_avgs.Properties.VariableNames =   {'GeneID', 'x_0_25', 'x0_25', 'x1', 'x3', 'x6', 'x9', 'x11', 'x11_75', 'x12_25', 'x13', 'x18', 'x23'};
writetable(nonlinear_rpkm_avgs, 'nonlinear_rpkm_avgs.xlsx')