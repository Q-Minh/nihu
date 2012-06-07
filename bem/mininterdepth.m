function mindepth = mininterdepth(tree, symm)
%MININTERDEPTH  Find minimal depth with nonempty interaction list
% MINDEPTH = MININTERDEPTH(TREE, SYMM) returns the minimal depth index in a
% cluster tree where the interaction list or the image interaction list is
% not empty

%% Argument check and default parameters
% narginchk(1,2);

if nargin == 0
    symm = 0;
end

%%
depth = length(tree)-1; % tree depth
mindepth = depth + 2;
for iL = 1 : depth + 1
    if ~isempty(tree(iL).interlist)
        mindepth = iL;
        break;
    end
    if symm
        if ~isempty(tree(iL).iminterlist)
            mindepth = iL;
            break;
        end
    end
end

end