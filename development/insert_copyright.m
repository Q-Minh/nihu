function insert_copyright(filein)
%INSERT_COPYRIGHT Insert Copyright rows in matlab m file
%  INSERT_COPYRIGHT(FILENAME) inserts copyright rows in the m file
%  FILENAME.

%   Copyright 2008-2010 P. Fiala
%   Budapest University of Technology and Economics
%   Dept. of Telecommunications

%% Argument check
error(nargchk(0,1,nargin));

%% 
if nargin == 0  % call function for all m files in directory
    d = dir('*.m');
    for i = 1 : length(d)
        insert_copyright(d(i).name);
    end
else
    %% Initialzation
    texttowrite  = {'%   Copyright 2008-2010 P. Fiala'
        '%   Budapest University of Technology and Economics'
        '%   Dept. of Telecommunications'};
    texttoremove  = {'%   Copyright 2008-2010 P. Fiala'
        '%   Budapest University of Technology and Economics'
        '%   Dept. of Telecommunications'
        '% Peter Fiala'
        '% 2009'};
    state = 's';
    write = false;
    copy = true;

    %% Open m file and temporary output file
    fin = fopen(filein, 'rt');
    if fin == -1
        error('NiHu:insert_copyright:fopen',...
            'Could not open m-file %s', filein);
    end
    tempname = 'temp.tmp';
    fout = fopen(tempname, 'wt');
    if fout == -1
        fclose(fin);
        error('NiHu:insert_copyright:fopen',...
            'Could not open temporary file to write');
    end
    %% Copy content from m-file to output file
    while true
        line = fgetl(fin);  % read line and exit if end of file
        if line == -1
            break;
        end

        switch state
            case 's'    % start state
                if strncmp(line, '%', 1)
                    % if comment line, go to comment state
                    state = 'c'; %#ok<*NASGU>
                elseif ~strncmp(strtrim(line), 'function', 8)
                    % if not comment and not 'function' line, write
                    % copyright content and go to text state
                    write = true;
                    state = 't';
                end
            case 'c'    % comment state
                if ~strncmp(line, '%', 1)
                    % if not comment line, print copyright and go to text state
                    write = true;
                    state = 't';
                end
            case 't'    % text state
                if copy || ~isempty(line)
                    for i = 1 : length(texttoremove)
                        if strncmp(line, texttoremove{i}, length(texttoremove{i}))
                            copy = false;
                            break;
                        end
                        copy = true;
                    end
                end
        end
        
        % enter copyright if requested
        if write
            fprintf(fout, '\n');
            for i = 1 : length(texttowrite)
                fprintf(fout, '%s\n', texttowrite{i});
            end
            write = false;
        end
        
        % copy file content if requested
        if copy
            fprintf(fout, '%s\n', line);
        end
    end
    
    %% Close input and output files and rename temp file to input
    fclose(fin);
    fclose(fout);
    if ~movefile(tempname, filein)
        error('NiHu:insert_copyright:fileWrite', ...
            'Cannot overwite file %s', filein);
    end
end

end
