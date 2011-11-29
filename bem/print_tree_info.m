function print_tree_info(tree)

disp(sprintf('Tree consists of %d levels (level 0 to %d)' , length(tree), length(tree)-1));
for l = 0 : length(tree)-1
    T = tree(l+1);
    disp(sprintf('Level %d:', l));
    disp(sprintf('\tCluster diameter: %f', T.diameter));
    disp(sprintf('\tNumber of Clusters: %d', size(T.coord,1)));
    if l > 0
        disp(sprintf('\tMax. length of nearfield list: %d', size(T.nearfield,2)));
    end
    if l > 1
        disp(sprintf('\tMax. length of interaction list: %d', size(T.interlist,2)));
        disp(sprintf('\tNumber of unique transfer distances: %d', size(T.D,1)));
    end
    if l == length(tree)-1
        disp(sprintf('\tMax. number of source nodes in leaf cluster: %d', size(T.nodsou,2)));
        disp(sprintf('\tMax. number of receiver nodes in leaf cluster: %d', size(T.nodrec,2)));
    end
end