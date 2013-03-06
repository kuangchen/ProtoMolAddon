#include "Hermite.h"

#include <protomol/base/MathUtilities.h>
#include <protomol/base/Report.h>
using namespace ProtoMol::Report;
using std::string;

namespace ProtoMol {
  //#define DEBUG_HERMITE

  //_________________________________________________________________ Hermite

  const string Hermite::keyword("Hermite");

  Hermite::Hermite():myInterOrder(0),
		     theta(NULL),
		     dTheta(NULL){
  }

  Hermite::Hermite(unsigned int order):myInterOrder(order),
				       theta(new Real[order]),
				       dTheta(new Real[order]){
  }

  Hermite::Hermite(unsigned int order, Real w):myInterOrder(order),
					       theta(new Real[order]),
					       dTheta(new Real[order]){
    set(w);
  }

  Hermite::~Hermite(){
    if(theta != NULL){ 
      delete [] theta;
      delete [] dTheta;
    }
  }

  Hermite::Hermite(const Hermite& Hermite){
    myInterOrder = Hermite.myInterOrder;
    theta  = new Real[myInterOrder];
    dTheta = new Real[myInterOrder];
    for(unsigned int k=0;k<myInterOrder;k++){
      theta[k] = Hermite.theta[k];
      dTheta[k] = Hermite.dTheta[k];      
    }
  }

  void Hermite::setOrder(unsigned int order){
    if(order == myInterOrder && theta != NULL) 
      return;
    delete [] theta;
    delete [] dTheta;
    theta = new Real[order];
    dTheta = new Real[order];
    myInterOrder = order;
  }

  void Hermite::set(Real w){
    switch(myInterOrder){
    case 0:
      report << error << "[Hermite::set] Interpolation order is zero!"<<endr;
      break;
    case 4:
      theta[0]  = (0.5*(1-w)*(1-w) * ((1-w)-1));     
      theta[1]  = (w*w*(1.5*w - 2.5) + 1);      
      theta[2]  = ((1-w)*(1-w)*(1.5*(1-w) - 2.5) + 1);
      theta[3]  = (0.5*w*w * (w-1));         
      dTheta[0] = (-(1.5*(1-w) - 1) * (1-w));           
      dTheta[1] = ((4.5*w - 5) * w);          
      dTheta[2] = (-((4.5*(1-w) - 5) * (1-w)));        
      dTheta[3] = (-(-(1.5*w - 1) * w));
      break;
    case 6:
      theta[0]  = (-1.0/24.0 * (1-(-w)) * (-w) * (1+(-w)) * (1+(-w)) * ((-w)+2));
      theta[1]  = (1.0/24.0 * (-w) *(1+(-w)) * (2+(-w)) * ((-5.0*(-w)+1)*(-w)+8));
      theta[2]  = (w*w-1) * (w-2) * ((-5.0/12.0*w+0.25)*w+0.5);
      theta[3]  = ((1-w)*(1-w)-1) * ((1-w)-2) * ((-5.0/12.0*(1-w)+0.25)*(1-w)+0.5);
      theta[4]  = (1.0/24.0 * (w-1) *(1+(w-1)) * (2+(w-1)) * ((-5.0*(w-1)+1)*(w-1)+8));
      theta[5]  = (-1.0/24.0 * (1-(w-1)) * (w-1) * (1+(w-1)) * (1+(w-1)) * ((w-1)+2));
      dTheta[0] = ((-1.0/24.0) * ((((5.0*(-w)+12)*(-w)+3)*(-w)-6)*(-w)-2));
      dTheta[1] = ((1.0/24.0) * ((((25.0*(-w)+56)*(-w)-3)*(-w)-52)*(-w)-16));
      dTheta[2] = (-(w*1.0/12.0) * (((25.0*w- 52)*w-15) * w+50));
      dTheta[3] = (((1-w)*1.0/12.0) * (((25.0*(1-w)- 52)*(1-w)-15) * (1-w)+50));
      dTheta[4] = (-(1.0/24.0) * ((((25.0*(w-1)+56)*(w-1)-3)*(w-1)-52)*(w-1)-16));
      dTheta[5] = (-(-1.0/24.0) * ((((5.0*(w-1)+12)*(w-1)+3)*(w-1)-6)*(w-1)-2));
      break;
    default:
      report << error << "[Hermite::set] Order "<<myInterOrder<<" not supported. Order 4 (qubic) and 6 (quintic) supported."<<endr;
      break;
    }
#if defined(DEBUG_HERMITE)
    report <<plain << "Hermite: order="<<myInterOrder<<", w="<<w<<endr;
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
