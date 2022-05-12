#ifndef __SMART_ASSERT
#define __SMART_ASSERT

//#define NDEBUG


void __M_Assert(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);
void __M_Error(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);
void __M_Warning(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);

#ifndef NDEBUG
#   define M_Assert_Break(Expr, Msg, breakExpr)\
    ({\
        if (Expr) {\
        __M_Error((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
        breakExpr;\
        }\
    })

#   define M_Assert_BreakSaveCheck(Expr, Msg, breakExpr)\
    ({\
        if (Expr) {\
        __M_Error((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
        breakExpr;\
        }\
    })

#   define M_Assert_Warning(Expr, Msg)\
    ({\
        if (Expr) {\
        __M_Warning((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
        }\
    })
#else
#   define M_Assert_Break(Expr, Msg, breakExpr);

#   define M_Assert_BreakSaveCheck(Expr, Msg, breakExpr)\
    ({\
        if(Expr) {\
        breakExpr;\
        }\
    })

#   define M_Assert_Warning(Expr, Msg);
#endif


#define PRINT_ONCE(...) ({ \
                        static int __printed = 0; \
                        if(!__printed) { \
                            printf(__VA_ARGS__); \
                            __printed=1; \
                        } \
                        })

#ifndef NULL
    #define NULL ((void *)0)
#endif

#endif // __SMART_ASSERT
