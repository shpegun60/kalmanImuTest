#include "smart_assert.h"
#include <assert.h>
#include <stdio.h>

void __M_Assert(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg)
{
    if (expr) {
        fprintf(stderr, "Assert failed:\t %s \nExpression:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
        //printf("Assert failed:\t %s:\nExpected:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
        abort();
    }
}

void __M_BreakMsg(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg)
{
    (void)expr;
    fprintf(stderr, "Assert failed:\t %s \nExpression:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
    //printf("Assert failed:\t %s:\nExpected:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
    abort();
}

void __M_Msg(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg)
{
    (void)expr;
    fprintf(stderr, "Assert failed:\t %s \nExpression:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
    //printf("Assert failed:\t %s:\nExpected:\t%s\nSource:\t\t%s, line %d\n", msg, expr_str,file, line);
}

