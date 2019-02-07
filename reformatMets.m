%% reformat metabolites in model
% Jeremy Primus
% February 5, 2019

%% description
% script to replace compartment designations (ex. [c]) with BiGG
% formatted identifiers (ex. _c)

function new_format_model = reformatMets(model)
    old_mets = model.mets;
    comp_index = cellfun(@(x) strfind(x, '['), old_mets);
    for i=1:length(old_mets)
        new_mets{i} = strcat(old_mets{i}(1:comp_index(i)-1), '_', old_mets{i}(comp_index(i)+1));
    end
    new_format_model = model;
    new_format_model.mets = new_mets';
end