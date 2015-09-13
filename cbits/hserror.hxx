#ifndef HS_ERROR_HXX
#define HS_ERROR_HXX

#include "HsFFI.h"

#ifdef __cplusplus
namespace HsMultivac {

  class Error
  {
  public:
    Error(HsStablePtr e) : m_err(e) {};
    HsStablePtr GetError () const
    {
      return m_err;
    }
    virtual ~Error() {};
  private:
    HsStablePtr m_err;
  };
}

#define CHECK_HS_ERROR(F) {\
  HsStablePtr err = F;\
  if (err) {\
    throw HsMultivac::Error(err);\
  }\
}

#endif

#endif
