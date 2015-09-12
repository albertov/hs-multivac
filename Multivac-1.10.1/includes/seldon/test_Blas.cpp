#define SELDON_WITH_CBLAS
#include "Seldon.hxx"
using namespace Seldon;

#include <iostream>
using namespace std;

int main()
{

  cout << "Seldon: compilation test with Blas" << endl;

  Vector<double> U(3), V(3);
  U.Fill(1.3);
  V.Fill();

  cout << "U = " << U << endl;
  cout << "V = " << V << endl;

  // 2.0 * U + V -> V
  Add(2.0, U, V);

  // Note: 'Add' calls the Blas function 'daxpy'.

  cout << "U + V = " << V << endl;

  return 0;

}
