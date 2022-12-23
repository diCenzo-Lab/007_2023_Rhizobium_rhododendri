function merge_rho62_1078

% Load data
temp3 = table2cell(readtable('temp3.txt'));
temp4 = table2cell(readtable('temp4.txt'));
orthologs = table2cell(readtable('orthologs.rho62.1078.txt', 'delimiter', ','));

% Merge the files
output = {};
for n = 1:length(orthologs)
    n
    pos = strmatch(orthologs{n,1}, temp3(:,1));
    pos2 = strmatch(orthologs{n,2}, temp4(:,1));
    output{n,1} = temp3{pos,2};
    output{n,2} = temp3{pos,3};
    output{n,3} = temp3{pos,4};
    temp = strsplit(temp4{pos2,2}, '(');
    output{n,4} = temp4{pos2,2};
    if strmatch('1078_Chromosome', temp)
        output{n,5} = 3664408 - temp4{pos2,3} + 1;
        output{n,6} = 3664408 - temp4{pos2,4} + 1;
    elseif strmatch('1078_ec1', temp)
        output{n,5} = 834411 - temp4{pos2,3} + 1;
        output{n,6} = 834411 - temp4{pos2,4} + 1;
    elseif strmatch('1078_ec2_pTi', temp)
        output{n,5} = 439071 - temp4{pos2,3} + 1;
        output{n,6} = 439071 - temp4{pos2,4} + 1;
    elseif strmatch('1078_ec3', temp)
        output{n,5} = 432998 - temp4{pos2,3} + 1;
        output{n,6} = 432998 - temp4{pos2,4} + 1;
    elseif strmatch('1078_ec4_Chromid', temp)
        output{n,5} = 304572 - temp4{pos2,3} + 1;
        output{n,6} = 304572 - temp4{pos2,4} + 1;
    elseif strmatch('1078_ec5', temp)
        output{n,5} = 302267 - temp4{pos2,3} + 1;
        output{n,6} = 302267 - temp4{pos2,4} + 1;
    end
    if strmatch('1078_Chromosome', output(n,4))
        output{n,7} = 'color=lgrey';
    elseif strmatch('1078_ec1', output(n,4))
        output{n,7} = 'color=vvdyellow';
    elseif strmatch('1078_ec2_pTi', output(n,4))
        output{n,7} = 'color=vdred';
    elseif strmatch('1078_ec3', output(n,4))
        output{n,7} = 'color=vdblue';
    elseif strmatch('1078_ec4_Chromid', output(n,4))
        output{n,7} = 'color=vdpurple';
    elseif strmatch('1078_ec5', output(n,4))
        output{n,7} = 'color=dorange';
    end
end

% Save the file
output_table = cell2table(output);
writetable(output_table, 'temp5.txt');
exit

