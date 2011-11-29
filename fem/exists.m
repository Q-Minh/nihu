function bool = exists(struct, fieldname)

if isfield(struct, fieldname)
    bool = ~isempty(getfield(struct, fieldname));
else
    bool = 0;
end