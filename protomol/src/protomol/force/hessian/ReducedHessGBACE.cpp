#include <protomol/force/hessian/ReducedHessGBACE.h>
#include <protomol/topology/GenericTopology.h>

#include <protomol/base/Report.h>

#include <iomanip>
#include <iostream>

#define PI 3.14169

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace std;

Matrix3By3 ReducedHessGBACE::operator()(Real a,
                                          const Vector3D &rij,
                                          const GenericTopology *topo,
                                          int atom1, int atom2,
                                          Real solvationparam, Real watersphereradius,
                                          ExclusionClass excl) const {
   
   Real dist = sqrt(a);

   //derivatives of burial term
   //Real btderv_ij = topo->atoms[atom1].myGBSA_T->btDerv1[atom2];
   //Real btderv_ji = topo->atoms[atom2].myGBSA_T->btDerv1[atom1];

   //cout <<setprecision(10) <<"HessGBACE :Atom "<<atom1<<", Atom2 "<<atom2<<", btderv_ij "<<btderv_ij<<", btderv_ji "<<btderv_ji<<endl;

   //Born radius
   Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
   Real bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;

   //derivatives of born radius
   //Real bornRadiusDerivative_ij = topo->atoms[atom1].myGBSA_T->bornRadiusDerivatives[atom2];
   //Real bornRadiusDerivative_ji = topo->atoms[atom2].myGBSA_T->bornRadiusDerivatives[atom1];

   //cout << setprecision(10) <<"HessGBACE : Atom "<<atom1<<", Atom2 "<<atom2<<", bornRadiusDerivative_ij "<<bornRadiusDerivative_ij<<", bornRadiusDerivative_ji "<<bornRadiusDerivative_ji<<endl;

   Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
   Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

   //offset radii ({\tilde{\rho}_{j}})
   Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
   Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;

   //first derivative of \Psi_i
   //Real psiderv_i_ij = offsetRadius_i*btderv_ij;

   Real psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;

   Real tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;

   Real tanh_i = tanh(tanhparam_i);

   Real tanhparam_i_derv = topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i;

   //data required for the second derivative of the burial term
   Real Lij = topo->atoms[atom1].myGBSA_T->Lvalues[atom2];
   Real Uij = topo->atoms[atom1].myGBSA_T->Uvalues[atom2];

   Real invLij = 1/Lij;
   Real invUij = 1/Uij;

   //Derivatives for calculation of the derivative of the born radii
   Real dLijdrij, dUijdrij, dCijdrij;

     //Scaling factors
     Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;


   //Check Equations (11-13)
   if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
   else dLijdrij = 0;

   if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
   else dUijdrij = 0;

   if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*(1/Lij)*(1/Lij)*dLijdrij;
   else dCijdrij = 0;

    //tanhk derivative
     Real tanhparam_derv_i = (topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i);

     Real S_i_term = (S_i*offsetRadius_i)/dist;
     Real S_j_term = (S_j*offsetRadius_j)/dist;


   //First Derivative of burial term
   Real btderv_ij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;

   Real psiderv_i_ij = offsetRadius_i*btderv_ij;

   Real bornRadiusDerivative_ij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*btderv_ij;

   //second derivative of the burial term

   Real d2BTijdrij2 = power(invLij,3)*power(dLijdrij,2) - power(invUij,3)*power(dUijdrij,2) + 0.5*power(invLij,3)*dLijdrij - 0.5*power(invUij,3)*dUijdrij + (dist/4)*(3*power(invUij,4)*power(dUijdrij,2) - 3*power(invLij,4)*power(dLijdrij,2)) + 0.5*(1/power(dist,3))*log(Lij/Uij) + 0.5*invLij*(1/(dist*dist))*(Lij*invUij*dUijdrij - dLijdrij) + 0.5*(1/dist)*invLij*invUij*(Lij*invUij*dUijdrij*dUijdrij - dUijdrij*dLijdrij) + 0.25*(1/dist)*invLij*invLij*dLijdrij*(Lij*invUij*dUijdrij - dLijdrij) - 0.25*(1/dist)*invLij*invUij*dUijdrij*(Lij*invUij*dUijdrij - dLijdrij) + ((power(S_j,2)*power(offsetRadius_j,2))/(4*dist*dist*dist))*(invLij*invLij - invUij*invUij) - ((power(S_j,2)*power(offsetRadius_j,2))/(2*dist*dist))*(invUij*invUij*invUij*dUijdrij - invLij*invLij*invLij*dLijdrij) + ((power(S_j,2)*power(offsetRadius_j,2))/(4*dist))*(3*power(invLij,4)*power(dLijdrij,2) - 3*power(invUij,4)*power(dUijdrij,2));


  Real d2Psi_i_drij2 = d2BTijdrij2*offsetRadius_i;

  Real alpha = topo->alphaObc;
  Real beta = topo->betaObc;
  Real gamma = topo->gammaObc;

  //second derivative of born radius
  Real d2Ridrij2 =  2*(1 - tanh_i*tanh_i)*(1-tanh_i*tanh_i)*power(psiderv_i_ij,2)*power(tanhparam_i_derv,2)*(power(bornRad_i,3)/power(radius_i,2)) -2*power(bornRad_i,2)*tanh_i*(1 - tanh_i*tanh_i)*power(psiderv_i_ij,2)*power(tanhparam_i_derv,2)*(1/radius_i) + (power(bornRad_i,2)/radius_i)*(1 - power(tanh_i,2))*(alpha*d2Psi_i_drij2 - 2*beta*power(psiderv_i_ij,2)- 2*beta*psi_i*d2Psi_i_drij2) + (power(bornRad_i,2)/radius_i)*(1 - power(tanh_i,2))*(6*gamma*psi_i*power(psiderv_i_ij,2) + 3*gamma*power(psi_i,2)*d2Psi_i_drij2);

   //data required for the second derivative of the burial term
   Real Lji = topo->atoms[atom2].myGBSA_T->Lvalues[atom1];
   Real Uji = topo->atoms[atom2].myGBSA_T->Uvalues[atom1];

   Real invLji = 1/Lji;
   Real invUji = 1/Uji;


     //Derivatives for calculation of the derivative of the born radii
      Real dLjidrij, dUjidrij, dCjidrij;

      if (offsetRadius_j <= (dist - S_i*offsetRadius_i)) dLjidrij = 1;
      else dLjidrij = 0;

      if (offsetRadius_j < (dist + S_i*offsetRadius_i)) dUjidrij = 1;
      else dUjidrij = 0;

      if (offsetRadius_j <= (S_i*offsetRadius_i - dist)) dCjidrij = 2*(1/Lji)*(1/Lji)*dLjidrij;
      else dCjidrij = 0;

   Real psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;   
   Real tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;
   Real tanh_j = tanh(tanhparam_j);
   Real tanhparam_j_derv = topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j;

     //tanhi derivative
     Real tanhparam_derv_j = (topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j);

   //first derivative of burial term
   Real btderv_ji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;

   Real psiderv_j_ij = offsetRadius_j*btderv_ji;

   Real bornRadiusDerivative_ji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*btderv_ji;

   //second derivative of the burial term

   Real d2BTjidrji2 = power(invLji,3)*power(dLjidrij,2) - power(invUji,3)*power(dUjidrij,2) + 0.5*power(invLji,3)*dLjidrij - 0.5*power(invUji,3)*dUjidrij + (dist/4)*(3*power(invUji,4)*power(dUjidrij,2) - 3*power(invLji,4)*power(dLjidrij,2)) + 0.5*(1/power(dist,3))*log(Lji/Uji) + 0.5*invLji*(1/(dist*dist))*(Lji*invUji*dUjidrij - dLjidrij) + 0.5*(1/dist)*invLji*invUji*(Lji*invUji*dUjidrij*dUjidrij - dUjidrij*dLjidrij) + 0.25*(1/dist)*invLji*invLji*dLjidrij*(Lji*invUji*dUjidrij - dLjidrij) - 0.25*(1/dist)*invLji*invUji*dUjidrij*
(Lji*invUji*dUjidrij - dLjidrij) + ((power(S_i,2)*power(offsetRadius_i,2))/(4*dist*dist*dist))*(invLji*invLji - invUji*invUji) - ((power(S_i,2)*power(offsetRadius_i,2))/(2*dist*dist))*(invUji*invUji*invUji*dUjidrij - invLji*invLji*invLji*dLjidrij) + ((power(S_i,2)*power(offsetRadius_i,2))/(4*dist))*(3*power(invLji,4)*power(dLjidrij,2) - 3*power(invUji,4)*power(dUjidrij,2));


   Real d2Psi_j_drij2 = d2BTjidrji2*offsetRadius_j;

  //second derivative of born radius
  Real d2Rjdrij2 =  2*(1 - tanh_j*tanh_j)*(1-tanh_j*tanh_j)*power(psiderv_j_ij,2)*power(tanhparam_j_derv,2)*(power(bornRad_j,3)/power(radius_j,2)) -2*power(bornRad_j,2)*tanh_j*(1 - tanh_j*tanh_j)*power(psiderv_j_ij,2)*power(tanhparam_j_derv,2)*(1/radius_j) + (power(bornRad_j,2)/radius_j)*(1 - power(tanh_j,2))*(alpha*d2Psi_j_drij2 - 2*beta*power(psiderv_j_ij,2)- 2*beta*psi_j*d2Psi_j_drij2) + (power(bornRad_j,2)/radius_j)*(1 - power(tanh_j,2))*(6*gamma*psi_j*power(psiderv_j_ij,2) + 3*gamma*power(psi_j,2)*d2Psi_j_drij2);

   //second derivative of the burial term
   Matrix3By3 vec_rij_ij(rij, rij);

   //Real first_term = (bornRadiusDerivative_ij + bornRadiusDerivative_ji)*(1/dist);
   Real rad_sum_i = radius_i + watersphereradius;
   Real rad_sum_j = radius_j + watersphereradius;

   Real first_term = -24*PI*solvationparam*(power(rad_sum_i,2)*power(radius_i,6)*(1/power(bornRad_i,7))*bornRadiusDerivative_ij + power(rad_sum_j,2)*power(radius_j,6)*(1/power(bornRad_j,7))*bornRadiusDerivative_ji)*(1/dist);

   Matrix3By3 I(1, 0 ,0 , 0 , 1, 0, 0 , 0, 1);

   Matrix3By3 H((I - vec_rij_ij*(1/(dist*dist)))*first_term);

   //Real second_term =  (d2Ridrij2 + d2Rjdrij2);
   Real second_term = 24*PI*solvationparam*(power(rad_sum_i,2)*power(radius_i,6)*(1/power(bornRad_i,7))*((7.0/bornRad_i)*bornRadiusDerivative_ij*bornRadiusDerivative_ij - d2Ridrij2) + power(rad_sum_j,2)*power(radius_j,6)*(1/power(bornRad_j,7))*((7.0/bornRad_j)*bornRadiusDerivative_ji*bornRadiusDerivative_ji - d2Rjdrij2));

   H += vec_rij_ij*(1/(dist*dist))*second_term;

   return H;
   
}
