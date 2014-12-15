function plot_cluster(CTree, Blocks, c)

h = plot_box(CTree(c).bb);

neighbours = [
    Blocks(Blocks(:,1)==c,2)
    Blocks(Blocks(:,2)==c,1)
    ];
for n = 1 : size(neighbours,1)
    h = plot_box(CTree(neighbours(n)).bb);
end

end

function h = plot_box(bb)

D = max(diff(bb,1));
x0 = mean(bb,1);

p = [
    1 2 3 4
    5 6 7 8
    1 2 6 5
    2 3 7 6
    3 4 8 7
    4 1 5 8
    ];

q = [
    -1 -1 -1
    1 -1 -1
    1 1 -1
    -1 1 -1
    -1 -1 1
    1 -1 1
    1 1 1
    -1 1 1
    ] * D/2;

x = q(:,1);
y = q(:,2);
z = q(:,3);
h = patch(x(p')+x0(1), y(p')+x0(2), z(p')+x0(3), ones(size(p')));

end
