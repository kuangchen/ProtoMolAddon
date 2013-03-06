#include "BSpline.h"

#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>
using namespace ProtoMol::Report;
using std::string;

//#define DEBUG_BSPLINE
namespace ProtoMol {
  //_________________________________________________________________ 

  const string BSpline::keyword("BSpline");

  BSpline::BSpline():myInterOrder(0),theta(NULL),dTheta(NULL){
  }

  BSpline::BSpline(unsigned int order):myInterOrder(order),
				       theta(new Real[order]),
				       dTheta(new Real[order]){
  }

  BSpline::BSpline(unsigned int order, Real w):myInterOrder(order),
					       theta(new Real[order]),
					       dTheta(new Real[order]){
    set(w);
  }

  BSpline::~BSpline(){
    if(theta != NULL){ 
      delete [] theta;
      delete [] dTheta;
    }
  }

  BSpline::BSpline(const BSpline& bspline){
    myInterOrder = bspline.myInterOrder;
    theta  = new Real[myInterOrder];
    dTheta = new Real[myInterOrder];
    for(unsigned int k=0;k<myInterOrder;k++){
      theta[k] = bspline.theta[k];
      dTheta[k] = bspline.dTheta[k];      
    }
  }

  void BSpline::setOrder(unsigned int order){
    if(order == myInterOrder && theta != NULL) 
      return;
    delete [] theta;
    delete [] dTheta;
    theta = new Real[order];
    dTheta = new Real[order];
    myInterOrder = order;
  }

  void BSpline::set(Real w){
    switch(myInterOrder){
    case 0:
      report << error << "[BSpline::set] Interpolation order is zero!"<<endr;
      break;
    case 1:
      theta[0]  = 1;
      dTheta[0] = 0;
      break;
    case 2:
      theta[1]  = w;
      theta[0]  = 1.0 - w;
      dTheta[1] = 1;
      dTheta[0] = -1;
      break;
    default:
      theta[myInterOrder-1] = 0.0;
      theta[1]     = w;
      theta[0]     = 1.0 - w;
      for(unsigned int k=2;k<myInterOrder-1;k++){
	Real div = 1.0 / k;
	theta[k] = div * w * theta[k-1];
	for(unsigned int j=1;j<k;j++)
	  theta[k-j] = div * ((w+j)*theta[k-j-1]+(k+1-j-w)*theta[k-j]);
	theta[0] = div * (1-w)*theta[0];
      }

      dTheta[0] = -theta[0];
      for(unsigned int j=1;j<=myInterOrder-1;j++)
	dTheta[j] = theta[j-1] - theta[j];
    
      Real div = 1.0 / (myInterOrder-1.0);
      theta[myInterOrder-1] = div * w * theta[myInterOrder-2];
      for(unsigned int j=1;j<myInterOrder-1;j++)
	theta[myInterOrder-1-j] = div * ((w+j)* theta[myInterOrder-j-2]+(myInterOrder-j-w)* theta[myInterOrder-1-j]);

      theta[0] = div * (1-w)* theta[0];
      break;
    }
#if defined(DEBUG_BSPLINE)
    report <<plain << "BSpline: order="<<myInterOrder<<", w="<<w<<endr;
    report << plain << "theta: ";
    for(unsigned int i=0;i<myInterOrder;i++)
      report << theta[i]<<" ";
    report << endr;
    report << plain << "dTheta: ";
    for(unsigned int i=0;i<myInterOrder;i++)
      report << dTheta[i]<<" ";
    report << endr;
#endif
  }
}
