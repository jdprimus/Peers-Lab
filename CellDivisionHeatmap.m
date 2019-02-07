%% Cell division genes heatmap
% Jeremy Primus
% February 5, 2019

%% Program information

% Description
% Beginnings of a script to create heatmaps of RNA-seq data
% modified from transcriptomefigures.m to test with new transcriptome data

% Dependancies
% - COBRA toolbox
% - Synechocystis genome-scale metabolic model provided by Ashok Prasad and
%   Chintan Joshi
% - findSubSystemOfGenes.m provided by Chintan Joshi
% - RNA-seq data: copied what appear to be the relevant data from "C:\Users\Peers Lab\Desktop\Jeremy - Spring 2019\Synechocystis\Graham_Leo\Graham_Leo\20170706_MasterAnnotationPCC6803.csv"

%% import data
clear all
load('cell_division_genes.mat')
rna_table = readtable("C:\Users\Peers Lab\Desktop\Jeremy - Spring 2019\Synechocystis\testRNAs.xlsx");
[goi, ia, ~] = intersect(rna_table.GeneID, cell_division_genes);

rna_data = table2array(rna_table(ia, 2:end));
rna_data(rna_data==0) = 0.1;        % remove zeros from data, as this will mess up fold change and log transformation
rna_data = rna_data./rna_data(:,1);                                                     % translate data to fold change
rna_log2_foldchange = log2(rna_data);                                                   % LOG2 transform

%% define custom colormap
[noOfGenes, ~] = size(rna_data);    % used to set resolution of color scale 

centerScale = 0;                    % value to center the color scale on 

maxColor = [1 0 0];                 % color for maximum data value: red
centerColor = [1 1 1];              % color for central data value: white
minColor = [0 0 1];                 % color for minimum data value: blue

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


    %% create figures
    % initialization
    xvalues = {'-0.25', '0.25', '1', '3', '6', '9', '11', '11.75', '12.25', '13', '18', '23'};
    yvalues = goi;
    lightbar = imread('editedlightcycle.png');              % image to overlay light conditions
    
    figure()
    % light condition subplot
    subplot('Position', [0.1 0.7 0.8 0.15])                
    xl = xlim;
    yl = ylim;
    imshow(lightbar, [], 'Xdata', 100*xl, 'Ydata', 10*yl)   % display image, stretched to fit subplot
    title('Cell Division Genes')
    % heatmap subplot
    subplot('Position', [0.1 0.3 0.8 0.4])                  
    h = heatmap(xvalues, yvalues, rna_log2_foldchange, 'CellLabelColor','none');  % display heatmap in subplot
    h.Colormap = mymap;
    colorbar
    caxis([-2.1, 3.4])                                      % scales each image by the 90% CI for the complete dataset
    xlabel('ZT Hours')
    ylabel('Locus Tag')

    % title(subsys_of_interest) % for manual selection of subsystem of interest