#ifndef HS_ERROR_HXX
#define HS_ERROR_HXX

#include "HsFFI.h"

typedef HsStablePtr HsException;

#ifdef __cplusplus
namespace HsMultivac {

  class Error
  {
  public:
    Error(HsException e) : m_err(e) {};
    HsException GetExceptionPtr () const
    {
      return m_err;
    }
    virtual ~Error() {};
  private:
    HsException m_err;
  };
}

#define CHECK_HS_ERROR(F) {\
  HsException err = F;\
  if (err) {\
    throw HsMultivac::Error(err);\
  }\
}

#define CATCH\
  }\
  catch (Seldon::Error& Err)\
  {\
    Err.What();\
    return multivacError ("Seldon error");\
  }\
  catch (Multivac::CError& Err)\
  {\
    Err.What();\
    return multivacError ("Multivac error");\
  }\
  catch (HsMultivac::Error& e) {\
    return e.GetExceptionPtr();\
  }\
  catch (std::exception& Err)\
  {\
    return multivacError (Err.what());\
  }\
  catch (std::string& str)\
  {\
    return multivacError (str.c_str());\
  }\
  catch (const char* str)\
  {\
    return multivacError (str);\
  }\
  catch(...)\
  {\
    return multivacError ("Unknown exception");\
  }\

#endif

#endif
