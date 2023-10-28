function partition = node_partitions(adj)

nNode = size(adj,1);

partition = zeros(nNode,1);
iPartition = 1;
base = find(partition==0, 1, 'first');

while ~isempty(base)
    while 1
        neighbours = find(any(adj(base,:),1));
        if length(neighbours) == length(base)
            break;
        end
        base = neighbours;
    end
    partition(base) = iPartition;
    iPartition = iPartition+1;
    base = find(partition==0, 1, 'first');
end
