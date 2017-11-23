function [path pathcell] = read_svgpath(fname)

str = parseXML(fname);

pathcell = {};

extract_paths(str);

actpos = [0 0];

path = zeros(0, 9);

for p = 1 : length(pathcell)
    first = true;
    split = strsplit(pathcell{p}, {' ', ','});
    while ~isempty(split)
        switch split{1}
            case 'M'
                k = 0;
                while true
                    actpos = [str2double(split{2*k+2}) str2double(split{2*k+3})];
                    if isnan(str2double(split(2*k+4)))
                        break;
                    end
                    k = k + 1;
                end
                split = split((2*k+4):end);
            case 'm'
                k = 0;
                while true
                    actpos = actpos + [str2double(split{2}) str2double(split{3})];
                    if 2*k+4 > length(split) || isnan(str2double(split(2*k+4)))
                        break;
                    end
                    k = k + 1;
                end
                split = split((2*k+4):end);
            case 'L'
                if (first)
                    start = actpos;
                    first = false;
                end
                dst = [str2double(split{2}) str2double(split{3})];
                path(end+1,1:5) = [1 actpos dst];
                actpos = dst;
                split = split(4:end);
            case 'H'
                if (first)
                    start = actpos;
                    first = false;
                end
                dst = [str2double(split{2}) actpos(2)];
                path(end+1,1:5) = [1 actpos dst];
                actpos = dst;
                split = split(3:end);
            case 'V'
                if (first)
                    start = actpos;
                    first = false;
                end
                dst = [actpos(1) str2double(split{2})];
                path(end+1,1:5) = [1 actpos dst];
                actpos = dst;
                split = split(3:end);
            case 'h'
                if (first)
                    start = actpos;
                    first = false;
                end
                dst = [actpos(1)+ str2double(split{2}),  actpos(2)];
                path(end+1,1:5) = [1 actpos dst];
                actpos = dst;
                split = split(3:end);
            case 'v'
                if (first)
                    start = actpos;
                    first = false;
                end
                dst = [actpos(1) actpos(2) + str2double(split{2})];
                path(end+1,1:5) = [1 actpos dst];
                actpos = dst;
                split = split(3:end);
            case 'C'
                if first
                    start = actpos;
                    first = false;
                end
                k = 0;
                while (true)
                    c1 = [str2double(split{6*k+2}) str2double(split{6*k+3})];
                    c2 = [str2double(split{6*k+4}) str2double(split{6*k+5})];
                    dst = [str2double(split{6*k+6}) str2double(split{6*k+7})];
                    path(end+1,1:9) = [2 actpos c1 c2 dst];
                    if isnan(str2double(split{6*k+8}))
                        break;
                    end
                    k = k+1;
                    actpos = dst;
                end
                split = split(6*k+8:end);
                
            case 'c'
                if first
                    start = actpos;
                    first = false;
                end
                k = 0;
                while (true)
                    dc1 = [str2double(split{6*k+2}) str2double(split{6*k+3})];
                    dc2 = [str2double(split{6*k+4}) str2double(split{6*k+5})];
                    ddst = [str2double(split{6*k+6}) str2double(split{6*k+7})];
                    c1 = actpos + dc1;
                    c2 = c1 + dc2;
                    dst = c2 + ddst;
                    path(end+1,1:9) = [2 actpos c1 c2 dst];
                    if (6*k+8 > length(split) || isnan(str2double(split{6*k+8})))
                        break;
                    end
                    k = k+1;
                    actpos = dst;
                end
                split = split(6*k+8:end);
            case {'Z', 'z'}
                path(end+1,1:5) = [1 actpos start];
                split = split(2:end);
            case ''
                split = split(2:end);
            otherwise
                error('Unprocessed svg path element %s', split{1});
        end
    end
end

    function extract_paths(str)
        for i = 1 : numel(str)
            if strcmp(str(i).Name, 'path')
                idx = strcmp({str.Attributes.Name}, 'd');
                pathcell{end+1} = str.Attributes(idx).Value;
            end
            
            for c = 1 : length(str(i).Children)
                extract_paths(str(i).Children(c));
            end
        end
    end
end % of function read_svg_path
