#ifndef MEMORY_SAVER_HXX


namespace HsMultivac
{


  template <class T>
  class CMemorySaver
  {
  public:

    CMemorySaver(vector<Curve<T> >& storage, T period)
      : m_fronts(storage)
      , m_period(period)
      , m_nextTime(0)
      {}

    void SaveAtTheBeginning(CMesh<T>& Mesh,
				    CSpeedFunction<T>& F,
				    CLevelSet<T>& Phi,
				    CInitializer<T>& Initializer) {};

    void SaveAtCurrentIteration(CMesh<T>& Mesh,
					CSpeedFunction<T>& F,
					CLevelSet<T>& Phi,
					T time, int iter,
					CInitializer<T>& Initializer)
    {
      if (time >= m_nextTime) {
        m_fronts.push_back(Curve<T>());
        Initializer.BuildCurveForDisplay(iter, Mesh, Phi);
        m_fronts.back().Copy(Initializer.GetFront());
        m_nextTime += m_period;
      }
    }

    void SaveAtTheEnd(CMesh<T>& Mesh,
			      CSpeedFunction<T>& F,
			      CLevelSet<T>& Phi,
			      Vector<T>& time, int iter,
			      CInitializer<T>& Initializer) {}

  private:
    vector<Curve<T> >& m_fronts;
    T                  m_period;
    T                  m_nextTime;
  };  // CMemorySaver.


}  // namespace HsMultivac.


#define MEMORY_SAVER_HXX
#endif
