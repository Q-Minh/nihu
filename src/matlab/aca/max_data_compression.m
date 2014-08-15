function comp = max_data_compression(rc, cc, bfar)

nTotal = dot(rc(bfar(:,1),2), cc(bfar(:,2),2), 1);
nMin = sum(2 * (rc(bfar(:,1),2) + rc(bfar(:,2),2)));
comp = nTotal/nMin;

end