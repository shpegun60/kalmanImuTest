#ifndef __SMART_ASSERT
#define __SMART_ASSERT

//#define NDEBUG


void __M_Assert(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);
void __M_BreakMsg(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);
void __M_Msg(const char* expr_str, unsigned char expr, const char* file, int line, const char* msg);

#ifndef NDEBUG
#   define M_Assert_Break(Expr, Msg, breakExpr) \
    if (Expr) {\
    __M_BreakMsg((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
     breakExpr;\
    }

#   define M_Assert_BreakSaveCheck(Expr, Msg, breakExpr) \
    if (Expr) {\
    __M_BreakMsg((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
    breakExpr;\
    }

#   define M_Assert_Warning(Expr, Msg, breakExpr) \
    if (Expr) {\
    __M_Msg((#Expr), (Expr), (__FILE__), (__LINE__), (Msg));\
    breakExpr;\
    }
#else
#   define M_Assert_Break(Expr, Msg, breakExpr);

#   define M_Assert_BreakSaveCheck(Expr, Msg, breakExpr) \
    if(Expr) {\
    breakExpr;\
    }

#   define M_Assert_Warning(Expr, Msg, breakExpr) \
    if(Expr) {\
    breakExpr;\
    }
#endif

#endif // __SMART_ASSERT
