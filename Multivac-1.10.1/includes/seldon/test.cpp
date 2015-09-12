#include "Seldon.hxx"
using namespace Seldon;

#include <iostream>
using namespace std;

int main()
{

  cout << "Seldon: compilation test" << endl;

  Vector<double> V(3);
  V.Fill();

  cout << "Vector: " << V << endl;

  return 0;

}
