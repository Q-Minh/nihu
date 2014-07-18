function democontents2xml(directory)
%DEMOCONTENTS2XML
%  SOME TEXT

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% open input and output files
fin = fopen(fullfile(directory, 'Contents.m'), 'r');
if fin == -1
    return;
end

fout = fopen(fullfile(directory,'demos.xml'), 'w');
if fout == -1
    fclose(fin);
    return;
end

%% print header
fprintf(fout, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fout, '<demos>');
fprintf(fout, '<name>NiHu</name>\n');
fprintf(fout, '<type>toolbox</type>\n');
fprintf(fout, '<icon>$toolbox/matlab/icons/matlabicon.gif</icon>\n');
fprintf(fout, '<description>These demos demonstrate how to use the toolbox NiHu.</description>\n');

insection = false;
while true
    line = fgetl(fin);
    if line == -1
        break;
    end
    if length(line) <= 2
        continue;
    end
    line = line(2:end); % remove comment sign
    % remove trailing white spaces
    while line(1) == ' '
        line = line(2:end);
    end
    if isempty(line)
        continue;
    end
    
    lim = strfind(line, ' - ');
    if isempty(lim)
        if insection
            fprintf(fout, '</demosection>\n');
            insection = false;
        end
        sectionname = line;
    else
        if ~insection
            fprintf(fout, '<demosection>\n<label>%s</label>', sectionname);
            insection = true;
        end
        mname = line(1:lim-1);
        while mname(end) == ' '
            mname = mname(1:end-1);
        end
        description = line(lim+3:end);
        fprintf(fout, '<demoitem>\n');
        fprintf(fout,'<label>%s</label>\n', description);
        fprintf(fout,'<callback>%s</callback>\n', mname);
        fprintf(fout,'<file>html/%s.html</file>\n', mname);
        fprintf(fout, '</demoitem>\n');
    end
end

if insection
    fprintf(fout, '</demosection>\n');
end

fprintf(fout, '</demos>');

fclose(fout);
fclose(fin);
end
