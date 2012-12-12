function answer = paramdlg(Title, Items)

nItems = size(Items,1);

answer = Items(:,2);

AnswerNeeded = 1;
while AnswerNeeded
    AnswerNeeded = 0;
    answer = inputdlg(Items(:,1), Title, 1, answer);
    if isempty(answer)
        AnswerNeeded = 0;
    else
        for iItems = 1 : nItems
            x = answer{iItems};
            ok = eval(Items{iItems,3});
            if ~ok
                errordlg(Items{iItems, 4}, 'Parameter Error', 'modal');
                uiwait;
                AnswerNeeded = 1;
                break;
            end
        end
    end
end