#include <iostream>
#include <fstream>
#include <LRSpline/LRSplineSurface.h>

using namespace std;
using namespace LR;

int main(int argc, char **argv) {
  // setup problem parameters
  int p1 = 3;   // polynomial order (degree + 1)
  int p2 = 3; 
  int n1 = 22;  // number of basisfunctions 
  int n2 = 22;  // remember that n=p gives a single element (knot span)
                // i.e. when changing p, also change n
      
  int lr_steps  = 3;    // Nr of refinement steps
  double radius = 0.05; // Radius for refinement around crack tip
        
  // construct object
  LRSplineSurface lr(n1,n2, p1,p2);

  // Refinement loop   
  for (int ref=0; ref<lr_steps; ++ref) {
    std::vector<int> idx;

    int i = 0;
    for(Basisfunction* b : lr.getAllBasisfunctions()) {
      // Find center of basis support
      double cx = (b->getParmin(0) + b->getParmax(0))/2.0;
      double cy = (b->getParmin(1) + b->getParmax(1))/2.0;

      // Refine if center is close to crack tip
      if ( (cx-0.5)*(cx-0.5) + (cy-0.5)*(cy-0.5) < radius*radius )
        idx.push_back(i);

      ++i;
    }

    lr.refineBasisFunction(idx);
  }

  // make the C^{-1} crack
  lr.insert_const_v_edge(0.5, // constant v-value
                         0,   // u_start
                         0.5, // u_end
                         p2); // multiplicity m (continuity is given by p-m-1)

  // write results to file
  ofstream myfile;
  myfile.open("unit_slit.lr", std::ofstream::trunc);
  myfile << lr << endl;
  myfile.close();
} 
