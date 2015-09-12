#ifndef SPEEDCALLBACK_CPP
#define SPEEDCALLBACK_CPP

#include "speedcallback.hpp"

namespace HsMultivac
{

  template <class T>
  CSpeedCallback<T>::CSpeedCallback(SpeedFuncCallbacks cb) throw()
    : m_cb(cb)
  {
    this->dependence_position = cb.dependencePosition;
    this->dependence_curvature = cb.dependenceCurvature;
    this->dependence_time = cb.dependenceTime;
    this->dependence_normal = cb.dependenceNormal;
  }

  template <class T>
  CSpeedCallback<T>::~CSpeedCallback()  throw() {}



  template <class T>
  inline void CSpeedCallback<T>::Init(CMesh<T>& Mesh)
  {
    this->Values.Reallocate(Mesh.GetNx(), Mesh.GetNy());
    this->Values.Fill(0.0);
  }


  template <class T>
  inline T CSpeedCallback<T>::operator() (T x, T y, T time) const
  {
    return (*m_cb.fastMarchSpeed)(x, y, time);
  }


  template <class T>
  inline T CSpeedCallback<T>::operator() (T x, T y, T time,
				      T nx, T ny, T curvature) const
  {
    return (*m_cb.narrowBandSpeed)(x, y, time, nx, ny, curvature);

  }

  template <class T>
  inline T CSpeedCallback<T>::GetMaxF1(T DxMin, T DxMax,
				   T DyMin, T DyMax,
				   T norm2) const
  {
    return (*m_cb.maxF1)(DxMin, DxMax, DyMin, DyMax, norm2);
  }


  template <class T>
  inline T CSpeedCallback<T>::GetMaxF2(T DxMin, T DxMax,
				   T DyMin, T DyMax,
				   T norm2) const
  {
    return (*m_cb.maxF2)(DxMin, DxMax, DyMin, DyMax, norm2);
  }

  template <class T>
  inline T CSpeedCallback<T>::GetDerivatives(T x, T y, T nx, T ny, T t,
					 T& dFdp, T& dFdx, T& dFdy,
					 T& dFdnx, T& dFdny) const
  {
    return 0;
  }

  template <class T>
  inline T CSpeedCallback<T>::Get2ndDerivatives(T x, T y, T nx, T ny, T t,
					    T& dFdpdp, T& dFdpdx, T& dFdpdy,
					    T& dFdpdnx, T& dFdpdny,
					    T& dFdxdx, T& dFdxdy,
					    T& dFdxdnx, T& dFdxdny,
					    T& dFdydy, T& dFdydnx,
					    T& dFdydny, T& dFdnxdnx,
					    T& dFdnxdny, T& dFdnydny) const
  {
    dFdpdp = 0.0;
    dFdpdx = 0.0;
    dFdpdy = 0.0;
    dFdpdnx = 0.0;
    dFdpdny = 0.0;
    dFdxdx = 0.0;
    dFdxdy = 0.0;
    dFdxdnx = 0.0;
    dFdxdny = 0.0;
    dFdydy = 0.0;
    dFdydnx = 0.0;
    dFdydny = 0.0;
    dFdnxdnx = 0.0;
    dFdnxdny = 0.0;
    dFdnydny = 0.0;
    return 0;
  }
}  // namespace Multivac.

#endif
