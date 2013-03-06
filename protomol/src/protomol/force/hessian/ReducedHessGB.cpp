#include <protomol/force/hessian/ReducedHessGB.h>
#include <protomol/topology/GenericTopology.h>

#include <protomol/base/Report.h>

#include <iomanip>
#include <iostream>

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace std;


Matrix3By3 ReducedHessGB::operator()( Real a,
                                          const Vector3D &rij,
                                          const GenericTopology *topo,
                                          int atom1, int atom2,
                                          int numatoms,
                                          Real soluteDielec, Real solventDielec,
                                          ExclusionClass excl) const {




   Real dist = sqrt(a);

   //derivatives of burial term


   //Born radius
   Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
   Real bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;



   //cout << setprecision(10) <<"HessGBACE : Atom "<<atom1<<", Atom2 "<<atom2<<", bornRadiusDerivative_ij "<<bornRadiusDerivative_ij<<", bornRadiusDerivative_ji "<<bornRadiusDerivative_ji<<endl;

   Real radius_i = topo->atoms[atom1].myGBSA_T->vanDerWaalRadius;
   Real radius_j = topo->atoms[atom2].myGBSA_T->vanDerWaalRadius;

   //offset radii ({\tilde{\rho}_{j}})
   Real offsetRadius_i = radius_i - topo->atoms[atom1].myGBSA_T->offsetRadius;
   Real offsetRadius_j = radius_j - topo->atoms[atom2].myGBSA_T->offsetRadius;

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


   //grab f_ij 
   //Real fGB_ij = topo->atoms[atom1].myGBSA_T->fij[atom2];
     Real expterm_ij = (dist*dist)/(4*bornRad_i*bornRad_j);

     Real exp_w = exp(-expterm_ij);
   Real fGB_ij = sqrt(dist*dist + bornRad_i*bornRad_j*exp_w);


   Real dfGBdrij = 0.5*(1/fGB_ij)*((2*dist - 0.5*exp(-expterm_ij)*dist) + exp(-expterm_ij)*bornRadiusDerivative_ij*(bornRad_j + (dist*dist)/(4*bornRad_i)) + exp(-expterm_ij)*bornRadiusDerivative_ji*(bornRad_i + (dist*dist)/(4*bornRad_j)));
    /*Real dfGBdrij = 0.5*(1/fGB_ij)*
    (
     2*dist + exp(-expterm_ij)*bornRadiusDerivative_ij*bornRad_j + exp(-expterm_ij)*bornRadiusDerivative_ji*bornRad_i
     + bornRad_i*bornRad_j*exp(-expterm_ij)*(
                                             -0.5*dist/(bornRad_i*bornRad_j) + 0.25*dist*dist/(bornRad_i*bornRad_i*bornRad_j)*bornRadiusDerivative_ij
                                             + 0.25*dist*dist/(bornRad_i*bornRad_j*bornRad_j)*bornRadiusDerivative_ji
                                             )
     );*/

     //second derivative
   Real d2fGBdrij2 = (-1/fGB_ij)*dfGBdrij*dfGBdrij + (1/fGB_ij) + 0.5*(1/fGB_ij)*exp_w*(d2Ridrij2*bornRad_j + 2*bornRadiusDerivative_ij*bornRadiusDerivative_ji + 2*bornRadiusDerivative_ij*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) + bornRad_i*d2Rjdrij2 + 2*bornRad_i*bornRadiusDerivative_ji*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji)+ bornRad_i*bornRad_j*(-(1/(2*bornRad_i*bornRad_j)) + ((dist)/(bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist)/(bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji - ((dist*dist)/(2*bornRad_i*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ij*bornRadiusDerivative_ji ) + bornRad_i*bornRad_j*( ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*d2Ridrij2 + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*d2Rjdrij2  - ((dist*dist)/(2*bornRad_i*power(bornRad_j,3)))*power(bornRadiusDerivative_ji,2) - ((dist*dist)/(2*bornRad_j*power(bornRad_i,3)))*power(bornRadiusDerivative_ij,2)) + bornRad_i*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji )*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) );

   //charges on atoms i and j
   Real charge_i = topo->atoms[atom1].scaledCharge;
   Real charge_j = topo->atoms[atom2].scaledCharge;

   Real second_derivative_born = charge_i*charge_j*( (1/power(fGB_ij,2))*d2fGBdrij2 - (2/power(fGB_ij,3))*dfGBdrij*dfGBdrij );

   // self terms
   second_derivative_born += power(charge_i, 2) / power(bornRad_i, 2) * (d2Ridrij2 / 2 - power( bornRadiusDerivative_ij, 2) / bornRad_i);
   second_derivative_born += power(charge_j, 2) / power(bornRad_j, 2) * (d2Rjdrij2 / 2 - power( bornRadiusDerivative_ji, 2) / bornRad_j);

   Real first_derivative_born = (charge_i*charge_j)*(1/power(fGB_ij,2))*dfGBdrij;

   //self terms
   first_derivative_born += power(charge_i, 2) * bornRadiusDerivative_ij / (2 * power(bornRad_i, 2));
   first_derivative_born += power(charge_j, 2) * bornRadiusDerivative_ji / (2 * power(bornRad_j, 2));

   Real d2Gikterm = 0;
   Real dGikterm = 0;
   Real d2Gjkterm = 0;
   Real dGjkterm = 0;
   //Real dfGBik, d2fGBik, dfGBjk, d2fGBjk;

  //summation of first derivatives
  //new N^2 handeling of interaction with other atoms
  if (!topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms) {
    topo->atoms[atom1].myGBSA_T->partialGBForceTerms = FirstDerivativeFGB(topo, atom1);
    topo->atoms[atom1].myGBSA_T->havePartialGBForceTerms = true;
  }

  dGikterm += ( topo->atoms[atom1].myGBSA_T->partialGBForceTerms - FirstDerivativeFGBSumError(topo, atom1, atom2) ) * bornRadiusDerivative_ij;
  
  if (!topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms) {
    topo->atoms[atom2].myGBSA_T->partialGBForceTerms = FirstDerivativeFGB(topo, atom2);
    topo->atoms[atom2].myGBSA_T->havePartialGBForceTerms = true;
  }

  dGjkterm += ( topo->atoms[atom2].myGBSA_T->partialGBForceTerms - FirstDerivativeFGBSumError(topo, atom2, atom1) )* bornRadiusDerivative_ji; 
  
  //summation of second derivatives
  //Atom i
  //new N^2 handeling of interaction with other atoms
  if (!topo->atoms[atom1].myGBSA_T->havePartialGBHessianTerms) {
    const secondDerivativeRawTerms sdrtatom = SecondDerivativeFGB(topo, atom1);
    topo->atoms[atom1].myGBSA_T->partialGBHessianTerms_term1 = sdrtatom.term1;
    topo->atoms[atom1].myGBSA_T->partialGBHessianTerms_term2 = sdrtatom.term2;
    topo->atoms[atom1].myGBSA_T->havePartialGBHessianTerms = true;
  }

  //const secondDerivativeRawTerms sdrtatom1 = SecondDerivativeFGB(topo, atom1);
  const secondDerivativeRawTerms sdrtatom1Err = SecondDerivativeFGBSumError(topo, atom1, atom2);
  d2Gikterm += (topo->atoms[atom1].myGBSA_T->partialGBHessianTerms_term1 - sdrtatom1Err.term1) * d2Ridrij2;
  d2Gikterm += (topo->atoms[atom1].myGBSA_T->partialGBHessianTerms_term2 - sdrtatom1Err.term2) * bornRadiusDerivative_ij * bornRadiusDerivative_ij;
  
  //Atom j
  //new N^2 handeling of interaction with other atoms
  if (!topo->atoms[atom2].myGBSA_T->havePartialGBHessianTerms) {
   const secondDerivativeRawTerms sdrtatom = SecondDerivativeFGB(topo, atom2);
   topo->atoms[atom2].myGBSA_T->partialGBHessianTerms_term1 = sdrtatom.term1;
   topo->atoms[atom2].myGBSA_T->partialGBHessianTerms_term2 = sdrtatom.term2;
   topo->atoms[atom2].myGBSA_T->havePartialGBHessianTerms = true;
  }
  
  //const secondDerivativeRawTerms sdrtatom2 = SecondDerivativeFGB(topo, atom2);
  const secondDerivativeRawTerms sdrtatom2Err = SecondDerivativeFGBSumError(topo, atom2, atom1);
  d2Gjkterm += (topo->atoms[atom2].myGBSA_T->partialGBHessianTerms_term1 - sdrtatom2Err.term1) * d2Rjdrij2;
  d2Gjkterm += (topo->atoms[atom2].myGBSA_T->partialGBHessianTerms_term2 - sdrtatom2Err.term2) * bornRadiusDerivative_ji * bornRadiusDerivative_ji;

  //finish 2nd derivative
   second_derivative_born += d2Gikterm + d2Gjkterm;
   second_derivative_born *= ((1/soluteDielec) - (1/solventDielec));

   first_derivative_born += dGikterm + dGjkterm;
   first_derivative_born *= ((1/soluteDielec) - (1/solventDielec))*(1/dist);
   
   Matrix3By3 I(1, 0 ,0 , 0 , 1, 0, 0 , 0, 1);
   Matrix3By3 vec_rij_ij(rij, rij);

   Matrix3By3 H((I - vec_rij_ij*(1/(dist*dist)))*first_derivative_born);

   H += vec_rij_ij*(1/(dist*dist))*second_derivative_born;

   return H;

}

Real ReducedHessGB::FirstDerivativeFGB(
                          const GenericTopology *topo,
                          const int atom1) const {

  //sum over derivatives
  Real dGBik = 0.0;

  //born radius and charge of i
  const Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  const Real charge_i = topo->atoms[atom1].scaledCharge;
   
  unsigned int sz = topo->atoms.size();

  for(unsigned k=0; k<sz; k++){
    
    if( k == (unsigned)atom1 ) continue;
      
    //grab r_{ik} / r_{jk}
    const Real dist = topo->atoms[atom1].myGBSA_T->distij[k];

    //born radius of k
    const Real bornRad_k = topo->atoms[k].myGBSA_T->bornRad;

    const Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
    const Real exp_ik = exp(-expfactor);
    Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);

    const Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*(bornRad_k + ((dist*dist)/(4*bornRad_i)));

    const Real charge_k = topo->atoms[k].scaledCharge;

    dGBik += (charge_i*charge_k)*(1/power(fGBik,2))*dfGBikdrij;
  }

  return dGBik;
}

Real ReducedHessGB::FirstDerivativeFGBSumError(
                                       const GenericTopology *topo,
                                       const int atom1, const int atom2) const {
  
  //born radius and charge of i
  const Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  const Real charge_i = topo->atoms[atom1].scaledCharge;
  
  //grab r_{ik} / r_{jk}
  const Real dist = topo->atoms[atom1].myGBSA_T->distij[atom2];
    
  //born radius of k
  const Real bornRad_k = topo->atoms[atom2].myGBSA_T->bornRad;
    
  const Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
  const Real exp_ik = exp(-expfactor);
  Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);
    
  const Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*(bornRad_k + ((dist*dist)/(4*bornRad_i)));
    
  const Real charge_k = topo->atoms[atom2].scaledCharge;
    
  Real dGBik = (charge_i*charge_k)*(1/power(fGBik,2))*dfGBikdrij;
  
  return dGBik;
}

ReducedHessGB::secondDerivativeRawTerms ReducedHessGB::SecondDerivativeFGB(
                          const GenericTopology *topo,
                          const int atom1) const {
  
  //sum of second derivatives
  secondDerivativeRawTerms sdrt;
  
  sdrt.term1 = sdrt.term2 = 0.0;
  
  //born radius and charge of i
  const Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  const Real charge_i = topo->atoms[atom1].scaledCharge;

  unsigned int sz = topo->atoms.size();
  
  for(unsigned k=0; k<sz; k++){
    
    if( k == (unsigned)atom1 ) continue; //|| k == atom2 

    //grab r_{ik} / r_{jk}
    const Real dist = topo->atoms[atom1].myGBSA_T->distij[k];
  
    //born radius of k
    const Real bornRad_k = topo->atoms[k].myGBSA_T->bornRad;
    
    //charge of k
    const Real charge_k = topo->atoms[k].scaledCharge;

    const Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
    const Real exp_ik = exp(-expfactor);
    const Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);

    //first derivative up to bornRadiusDerivative_ij
    const Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*(bornRad_k + ((dist*dist)/(4*bornRad_i))); //bornRadiusDerivative_ij*

    //second derivative up to bornRadiusDerivative_ij (terms 1 and 3) and d2Ridrij2 for term 1
    const Real d2fGBikdrij2_term1 = -(1/fGBik)*dfGBikdrij*dfGBikdrij;//    *bornRadiusDerivative_ij*bornRadiusDerivative_ij;
    const Real d2fGBikdrij2_term2 = (1/(2.0*fGBik))*exp_ik*(bornRad_k + ((dist*dist)/(4.0*bornRad_i))); //d2Ridrij2
    const Real d2fGBikdrij2_term3 = (1/(2.0*fGBik))*exp_ik*(power(dist,4)/(16*power(bornRad_i,3)*bornRad_k));// *bornRadiusDerivative_ij*bornRadiusDerivative_ij;


    sdrt.term1 += (charge_i*charge_k)*(1/power(fGBik,2))*d2fGBikdrij2_term2; //*d2Ridrij2;
    
    sdrt.term2 += (charge_i*charge_k)*( (-2/power(fGBik,3))*dfGBikdrij*dfGBikdrij
                + (1/power(fGBik,2))*( d2fGBikdrij2_term1 +d2fGBikdrij2_term3 )
                                   ); //*bornRadiusDerivative_ij*bornRadiusDerivative_ij;
    //
    
  }
   
  return sdrt;

}

ReducedHessGB::secondDerivativeRawTerms ReducedHessGB::SecondDerivativeFGBSumError(
                                                                           const GenericTopology *topo,
                                                                           const int atom1, const int atom2) const {
  
  //second derivatives
  secondDerivativeRawTerms sdrt;
    
  //born radius and charge of i
  const Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
  const Real charge_i = topo->atoms[atom1].scaledCharge;
  
  //grab r_{ik} / r_{jk}
  const Real dist = topo->atoms[atom1].myGBSA_T->distij[atom2];
    
  //born radius of k
  const Real bornRad_k = topo->atoms[atom2].myGBSA_T->bornRad;
    
  //charge of k
  const Real charge_k = topo->atoms[atom2].scaledCharge;
    
  const Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
  const Real exp_ik = exp(-expfactor);
  const Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);
    
  //first derivative up to bornRadiusDerivative_ij
  const Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*(bornRad_k + ((dist*dist)/(4*bornRad_i))); //bornRadiusDerivative_ij*
    
  //second derivative up to bornRadiusDerivative_ij (terms 1 and 3) and d2Ridrij2 for term 1
  const Real d2fGBikdrij2_term1 = -(1/fGBik)*dfGBikdrij*dfGBikdrij;//    *bornRadiusDerivative_ij*bornRadiusDerivative_ij;
  const Real d2fGBikdrij2_term2 = (1/(2.0*fGBik))*exp_ik*(bornRad_k + ((dist*dist)/(4.0*bornRad_i))); //d2Ridrij2
  const Real d2fGBikdrij2_term3 = (1/(2.0*fGBik))*exp_ik*(power(dist,4)/(16*power(bornRad_i,3)*bornRad_k));// *bornRadiusDerivative_ij*bornRadiusDerivative_ij;
    
    
  sdrt.term1 = (charge_i*charge_k)*(1/power(fGBik,2))*d2fGBikdrij2_term2; //*d2Ridrij2;
    
  sdrt.term2 = (charge_i*charge_k)*( (-2/power(fGBik,3))*dfGBikdrij*dfGBikdrij
                                       + (1/power(fGBik,2))*( d2fGBikdrij2_term1 +d2fGBikdrij2_term3 )
                                       ); //*bornRadiusDerivative_ij*bornRadiusDerivative_ij;
  //
  
  return sdrt;
  
}



