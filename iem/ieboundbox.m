function bb = ieboundbox(Map)
temp = Map(end/2+1:end,:)-Map(1:end/2,:);
hMap = Map(1:end/2,:);
temp(temp ~= 0) = Inf*temp(temp ~= 0);
temp(temp == 0) = hMap(temp == 0);
bb = [min([Map;temp]),max([Map;temp])];
end