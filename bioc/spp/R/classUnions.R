##definitions of class unions
setClassUnion("list_Or_NULL", c("list", "NULL"))
setClassUnion("list_Or_missing_Or_NULL", c("list", "missing", "NULL"))
setClassUnion("char_Or_NULL", c("character", "NULL"))
setClassUnion("char_Or_missing", c("character", "missing"))
setClassUnion("integer_Or_missing", c("integer", "missing"))
setClassUnion("numeric_Or_missing", c("numeric", "missing"))
setClassUnion("numeric_Or_NULL", c("numeric", "NULL"))
setClassUnion("numeric_Or_missing_Or_NULL", c("numeric", "missing", "NULL"))
setClassUnion("logical_Or_missing", c("logical", "missing"))
setClassUnion("df_Or_missing", c("data.frame", "missing"))
setClassUnion("df_Or_NULL", c("data.frame", "NULL"))
setClassUnion("function_Or_missing", c("function", "missing"))
##


