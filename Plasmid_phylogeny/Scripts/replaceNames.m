function replaceNames

% Load data
namesToReplace = table2cell(readtable('intermediaryFiles/names_to_replace.txt', 'delimiter', "\t", 'ReadVariableNames', false));
fasta = table2cell(readtable('intermediaryFiles/ParA_trimal.txt', 'delimiter', "\t", 'ReadVariableNames', false));

% Replace the names
for n = 1:length(namesToReplace)
    pos = strmatch(namesToReplace{n,2}, fasta(:), 'exact');
    if ~isempty(pos)
        fasta{pos} = namesToReplace{n,1};
    end
end

% Save and exit
fasta_new = cell2table(fasta);
writetable(fasta_new, 'intermediaryFiles/ParA_trimal_renamed.txt', 'WriteVariableNames', false);
exit
