#ifndef HS_ERROR_HPP
#define HS_ERROR_HPP

#include "HsFFI.h"

#ifdef __cplusplus
namespace HsMultivac {

  class Error
  {
  public:
    Error(HsStablePtr e) : m_err(e) {};
    virtual ~Error() {};
  private:
    HsStablePtr m_err;
  };
}

extern "C" {
#endif // __cplusplus

void MVThrowError(HsStablePtr e);

#ifdef __cplusplus
}
void MVThrowError(HsStablePtr e) {
  throw HsMultivac::Error(e);
}
#endif

#endif
