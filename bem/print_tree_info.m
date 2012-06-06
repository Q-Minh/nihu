function print_tree_info(tree)

fprintf(1, 'Tree consists of %d levels (level 0 to %d)\n' , length(tree), length(tree)-1);
for l = 0 : length(tree)-1
    T = tree(l+1);
    fprintf(1, 'Level %d:\n', l);
    fprintf(1, '\tCluster diameter: %f\n', T.diameter);
    fprintf(1, '\tNumber of Clusters: %d\n', size(T.coord,1));
    if l > 0
        fprintf(1, '\tMax. length of nearfield list: %d\n', size(T.nearfield,2));
    end
    if l > 1
        fprintf(1, '\tMax. length of interaction list: %d\n', size(T.interlist,2));
        fprintf(1, '\tNumber of unique transfer distances: %d\n', size(T.D,1));
    end
    if l == length(tree)-1
        fprintf(1, '\tMax. number of source nodes in leaf cluster: %d\n', size(T.nodsou,2));
        fprintf(1, '\tMax. number of receiver nodes in leaf cluster: %d\n', size(T.nodrec,2));
    end
end
