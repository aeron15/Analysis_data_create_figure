function idx=determine_index(list,query)

%DETERMINE_INDEX determines the index of a query in a cell array
% list is a cell array of strings
% query is a string

idx=find(strcmp(list,query));


end