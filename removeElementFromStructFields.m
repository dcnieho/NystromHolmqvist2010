function str = removeElementFromStructFields(str,idxOrBool)
str = structfun(@(x) removeElementFromVectors(x,idxOrBool),str, 'uni',false);


% helper func
function field = removeElementFromVectors(field,idxOrBool)

if ~isscalar(field)
    field(idxOrBool) = [];
end