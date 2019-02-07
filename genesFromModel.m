%% Genes from model
% Jeremy Primus
% January 16, 2019

%% Description
% code to extract gene sequences from model iJN678 at Bigg Models database
% planned updates: replace with Ashok's model, iSynCJ816
%% 
load('iJN678.mat')
url = 'http://bigg.ucsd.edu/models/iJN678/genes/';      % url for list of genes in model
for i = 1:length(iJN678.genes)                          % loop over genes
    geneurl = strcat(url, iJN678.genes{1});
    data = webread(geneurl);                            % read the webpage
    str = extractHTMLText(data);                        % this function is enabled by the Text Analytics Toolbox free trial
    gene_start = strfind(str, 'DNA Sequence') + length('DNA Sequence') + 1;
    gene_end = strfind(str, 'Protein Sequence') - 1;
    gene_seq{i} = strtrim(str(gene_start:gene_end));    % pull nucleotide sequence from webpage
end