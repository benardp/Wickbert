#include "trianglenearestpoint.h"
#include "distancetransform.h"
#include "utils.h"

int main(int argc, char **argv){

  SignedDistanceTransform::Test();

  TriangleNearestPoint::Test();
 
  TestUtils();

  return 0;
}
