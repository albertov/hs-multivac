#ifndef HS_ERROR_HPP
#define HS_ERROR_HPP

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

extern "C" {
#endif // __cplusplus

typedef void (*MVAction)();
void MVThrowError(HsStablePtr e);
HsStablePtr MVCatchError(MVAction);

#ifdef __cplusplus
}

void MVThrowError(HsStablePtr e) {
  throw HsMultivac::Error(e);
}

HsStablePtr MVCatchError(MVAction f)
{
  try {
    (*f)();
    return 0;
  } catch (HsMultivac::Error& e) {
    return e.GetError();
  }
}

#endif

#endif
