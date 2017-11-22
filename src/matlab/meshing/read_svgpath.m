function paths = read_svgpath(fname)

str = parseXML(fname);

paths = {};

extract_paths(str);

    function extract_paths(str)
        for i = 1 : numel(str)
            if strcmp(str(i).Name, 'path')
                idx = strcmp({str.Attributes.Name}, 'd');
                paths{end+1} = str.Attributes(idx).Value;
            end
            
            for c = 1 : length(str(i).Children)
                extract_paths(str(i).Children(c));
            end
        end
    end

end