#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Array.h>

using namespace std;
using namespace ProtoMol;
//____ Vector3DBlock


void Vector3DBlock::boundingbox(Vector3D &minbb, Vector3D &maxbb) const {
  if (size() <= 0) {
    minbb = Vector3D(Constant::MAXREAL, Constant::MAXREAL, Constant::MAXREAL);
    maxbb = Vector3D(-Constant::MAXREAL,
      -Constant::MAXREAL,
      -Constant::MAXREAL);
  }

  minbb = Vector3D(vec[0]);
  maxbb = Vector3D(vec[0]);

  const unsigned int count = size();
  for (unsigned int i = 1; i < count; ++i) {
    if (c[i*3] < minbb.c[0])
      minbb.c[0] = c[i*3];
    else if (c[i*3] > maxbb.c[0])
      maxbb.c[0] = c[i*3];
    if (c[i*3+1] < minbb.c[1])
      minbb.c[1] = c[i*3+1];
    else if (c[i*3+1] > maxbb.c[1])
      maxbb.c[1] = c[i*3+1];
    if (c[i*3+2] < minbb.c[2])
      minbb.c[2] = c[i*3+2];
    else if (c[i*3+2] > maxbb.c[2])
      maxbb.c[2] = c[i*3+2];
  }
}



Vector3D Vector3DBlock::sum() const {
  if (vec.empty())
    return Vector3D(0.0, 0.0, 0.0);

  //  Loop through all the atoms and remove the motion of center of mass
  //  using an advanced method -- Kahan's magic addition algorithm to get
  //  rid of round-off errors: Scientific Computing pp34.

  Vector3D sum(vec[0]);
  Vector3D tempC(0.0, 0.0, 0.0);
  for (unsigned int i = 1; i < size(); i++) {
    Vector3D tempX(vec[i]);
    Vector3D tempY(tempX - tempC);
    Vector3D tempT(sum + tempY);
    tempC = (tempT - sum) - tempY;
    sum = tempT;
  }

  return sum;
}

