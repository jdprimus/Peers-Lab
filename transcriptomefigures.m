%% Transcriptome data heatmap
% Jeremy Primus
% January 15, 2019

%% Description
% Beginnings of a script to create heatmaps of RNA-seq data

%% import data
rna_table = readtable("E:\RNAseq_PeersLab_EFRI\CLC_RNAseq_OutputDATA\20170413_RNAseqResultsAll_ReversMappingOnly\20170413_SummaryRNAseq12timpoints.xlsx");
rna_data = table2array(rna_table(1:50,3:10));
rna_data = rna_data./rna_data(:,1);     % translate data to fold change
rna_data(~isfinite(rna_data)) = 1;      % removes Inf values from data - should modify to treat initial reads of 0 differently
rna_log2_foldchange = log2(rna_data);   
%% define custom colormap
[noOfGenes, ~] = size(rna_data); 

centerScale = 0;        % value to center the color scale on 

maxColor = [1 0 0];     % color for maximum data value: red
centerColor = [1 1 1];  % color for central data value: white
minColor = [0 0 1];     % color for minimum data value: blue

% proportionally calculate where centerScale lies between min and max values
largest = max(max(rna_log2_foldchange));          
smallest = min(min(rna_log2_foldchange));
centerValue = noOfGenes*abs(centerScale - smallest)/(largest - smallest);

% Create color map ranging from centerColor to maxColor
mymap = [linspace(minColor(1),centerColor(1),centerValue)',...
         linspace(minColor(2),centerColor(2),centerValue)',...
         linspace(minColor(3),centerColor(3),centerValue)';...
         linspace(centerColor(1),maxColor(1),(noOfGenes-centerValue))',...
         linspace(centerColor(2),maxColor(2),(noOfGenes-centerValue))',...
         linspace(centerColor(3),maxColor(3),(noOfGenes-centerValue))'];

%% create heatmap
xvalues = rna_table.Properties.VariableNames(3:10);
yvalues = rna_table.Name(1:50);
h = heatmap(xvalues, yvalues, rna_log2_foldchange);
h.Colormap = mymap;