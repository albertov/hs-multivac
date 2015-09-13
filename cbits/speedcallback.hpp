#ifndef SPEEDCALLBACK_HPP
#define SPEEDCALLBACK_HPP

#include "hserror.hxx"

#ifdef __cplusplus

#include <cstdio>
#include "errors.cxx"

extern "C" {
#endif

/*
 * C interface
 */
typedef HsException (*FastMarchSpeedFunc)(double, double, double, double*);
typedef HsException (*NarrowBandSpeedFunc)(double, double, double, double, double, double, double*);
typedef HsException (*MaxFSpeedFunc)(double, double, double, double, double, double*);

typedef struct SpeedFuncCallbacks {
  FastMarchSpeedFunc fastMarchSpeed; 
  NarrowBandSpeedFunc narrowBandSpeed;
  MaxFSpeedFunc maxF1;
  MaxFSpeedFunc maxF2;
  int dependencePosition;
  int dependenceTime;
  int dependenceNormal;
  int dependenceCurvature;
} SpeedFuncCallbacks;

#ifdef __cplusplus
}


using namespace Multivac;

namespace HsMultivac
{

  

  template <class T>
  class CSpeedCallback: public CSpeedFunction<T>
  {


  public:

    CSpeedCallback(SpeedFuncCallbacks callbacks)  throw();
    ~CSpeedCallback()  throw();


    /***********
     * METHODS *
     ***********/

  public:
  
    virtual void Init(CMesh<T>& Mesh);

    virtual inline T operator() (T x, T y, T time) const;
    virtual inline T operator() (T x, T y, T time,
				 T nx, T ny, T curvature) const;

    virtual T GetMaxF1(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const;
    virtual T GetMaxF2(T Xmin, T Xmax, T Ymin, T Ymax, T norm2) const;

    virtual T GetDerivatives(T x, T y, T nx, T ny, T t,
			     T& dFdp, T& dFdx, T& dFdy,
			     T& dFdnx, T& dFdny) const;
    virtual T Get2ndDerivatives(T x, T y, T nx, T ny, T t,
				T& dFdpdp, T& dFdpdx, T& dFdpdy,
				T& dFdpdnx, T& dFdpdny,
				T& dFdxdx, T& dFdxdy,
				T& dFdxdnx, T& dFdxdny,
				T& dFdydy, T& dFdydnx,
				T& dFdydny, T& dFdnxdnx,
				T& dFdnxdny, T& dFdnydny) const;

  private:
    
    SpeedFuncCallbacks m_cb;

  };  // CSpeedCallback.


}  // namespace HsMultivac.


#endif // __cplusplus
#endif
