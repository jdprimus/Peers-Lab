%% Genes to Rxn Pathways
% Jeremy Primus
% January 17, 2019

%% Description
% code to map associated genes to reaction pathways and recreate figure
% from Ashok/Chintan's paper to verify the result 

%% Note: 1/22/2019
% this script produces results that differ from those cited in the paper.
% script 'findSubsystemOfGenes.m' provided by Chintan exactly reproduces
% the paper figure.  Not yet certain where the discrepancies are.
%%
clear all
j=1;
load('iCJ816_d.mat')                 % load model provided by Dr. Chintan Joshi
model2 = buildRxnGeneMat(iCJ816);    % matrix mapping genes to pathways
for i= 1:length(model2.rxnGeneMat(:,1))
    if i~=1 && any(strcmp(pathway, model2.subSystems{i}))               % if process already exists in pathway
        involved_genes{strcmp(pathway, model2.subSystems{i})} = vertcat(involved_genes{strcmp(pathway, model2.subSystems{i})}, model2.genes(model2.rxnGeneMat(i,:)==1));    % add genes to list of genes in respective pathway
    else    
        pathway{j} = model2.subSystems{i};                              % add new process to pathway
        involved_genes{j} = model2.genes(model2.rxnGeneMat(i,:)==1);    % add genes to new process  
        j = j+1;
    end
end

for j = 1:length(involved_genes)                                        % loop over list of lists
    involved_genes{j} = unique(involved_genes{j});                      % pull out the repeat genes in each list
    noOfinvolved_genes(j) = length(involved_genes{j});                  % determine the number of genes per pathway
end
[noOfinvolved_genes_sorted, index] = sort(noOfinvolved_genes);          % sort by the number of associated genes
pathway = categorical(pathway, pathway(index));                         % create pathways as categorical, order by number of associated genes
barh(pathway, noOfinvolved_genes)                                       % horizontal barplot of gene/pathway distribution