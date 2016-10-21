#include <iostream>
#include <fstream>
#include <LRSpline/LRSplineSurface.h>

using namespace std;
using namespace LR;

int main(int argc, char **argv) {
  // setup problem parameters
  int p1 = 2;   // polynomial order (degree + 1)
  int p2 = 2; 
  int n1 = 21;  // number of basisfunctions 
  int n2 = 21;  // remember that n=p gives a single element (knot span)
                // i.e. when changing p, also change n
  double lx = 4; // Length of domain
  double ly = 4;
      
  const int lr_steps = 3; // Nr of refinement steps
        
  // construct object
  LRSplineSurface lr = LRSplineSurface(n1,n2, p1,p2);

  for(Basisfunction* b : lr.getAllBasisfunctions()) {
    b->cp()[0] *= lx; // scale x-coordinate of this controlpoint
    b->cp()[1] *= ly; // scale y-coordinate of this controlpoint
  }

  if (lr_steps) {

    double delta = 0.14;
      
    for (int ref=0; ref<lr_steps; ++ref) {

      cout << "Refinement cycle #" << ref+1 << endl;

      std::vector<int> idx;

      int i = 0;
      for(Basisfunction* b : lr.getAllBasisfunctions()) {
	// Find center of basis support
	double cx = (b->getParmin(0) + b->getParmax(0))/2.0;
	double cy = (b->getParmin(1) + b->getParmax(1))/2.0;

	// Refine if center is close to line y=0.5 and x>0.5
	//if ( (cx > 0.45-delta && cx < 0.55+delta) && (cy < (0.5+delta) && cy > (0.5-delta)) ) {
	if ( cy < (0.5+delta) && cy > (0.5-delta) ) {
	  idx.push_back(i);
	  cout << "  Tag basis function " << i << " for refinement" << endl;
	}

	++i;
      }

      lr.refineBasisFunction(idx);

      delta = 0.67*delta;
	
    }

  }

  // write results to file
  ofstream myfile;
  myfile.open("sneddon.lr", std::ofstream::trunc);
  myfile << lr << endl;
  myfile.close();
} 
