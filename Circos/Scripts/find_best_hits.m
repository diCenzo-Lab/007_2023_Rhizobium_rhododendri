function find_best_hits

% Load data
rho62_1078 = table2cell(readtable('rho62_1078.txt'));
s1078_rho62 = table2cell(readtable('1078_rho62.txt'));
s932_1078 = table2cell(readtable('932_1078.txt'));
s1078_932 = table2cell(readtable('1078_932.txt'));

% Find blast-BBH for rho62 and 1078
first = rho62_1078;
second = s1078_rho62;
output_rho62_1078 = {};
for n = 1:length(first)
    n
    if(first{n,3} >= 50 && first{n,14} <= 1 * 10^-100)
        pos = strmatch(first(n,2), second(:,1), 'exact');
        pos2 = strmatch(second(pos,2), first(:,1), 'exact');
        if pos2 == n
            if(second{pos,3} >= 50 && second{pos,14} <= 1 * 10^-100)
                output_rho62_1078{end+1,1} = first{n,1};
                output_rho62_1078{end,2} = first{n,2};
            end
        end
    end
end

% Find blast-BBH for 932 and 1078
first = s932_1078;
second = s1078_932;
output_932_1078 = {};
for n = 1:length(first)
    n
    if(first{n,3} >= 50 && first{n,14} <= 1 * 10^-100)
        pos = strmatch(first(n,2), second(:,1), 'exact');
        pos2 = strmatch(second(pos,2), first(:,1), 'exact');
        if pos2 == n
            if(second{pos,3} >= 50 && second{pos,14} <= 1 * 10^-100)
                output_932_1078{end+1,1} = first{n,1};
                output_932_1078{end,2} = first{n,2};
            end
        end
    end
end

% Save
output_rho62_1078_table = cell2table(output_rho62_1078);
writetable(output_rho62_1078_table, 'orthologs.rho62.1078.txt');
output_932_1078_table = cell2table(output_932_1078);
writetable(output_932_1078_table, 'orthologs.932.1078.txt');
exit


