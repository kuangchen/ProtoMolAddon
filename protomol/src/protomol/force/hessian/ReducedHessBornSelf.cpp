#include <protomol/force/hessian/ReducedHessBornSelf.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/force/born/BornSwitch.h>

using namespace ProtoMol;
//____ ReducedHessBornSelf
// Calculates the Hessian of the Born self energy.
Matrix3By3 ReducedHessBornSelf::operator()(Real a,
                                          const Vector3D &rij,
                                          const GenericTopology *topo,
                                          int atom1, int atom2,
                                          int bornSwitch, Real dielecConst,
                                          ExclusionClass excl) const {

  Real dist = sqrt(a);
  Real rDist = 1.0 / dist;      
  // Atom 1 variables
  int type1 = topo->atoms[atom1].type;      // Atom variables
  int type2 = topo->atoms[atom2].type;
  //switches
  BornSwitch bSw(bornSwitch);  // Always quartic now  6/1/2008
  Real f_ij = bSw.switchValue(dist);
  Real df_ij = bSw.switchDerivative(dist);
  Real d2f_ij = bSw.switchSecondDerivative(dist);
  //Atom 1 Born Radius derivatives
  //**********************************
  Real C_i = topo->atomTypes[type1].mySCPISM_T->C_i;
  Real eta_i = topo->atoms[atom1].mySCPISM_A->eta;
  Real expCrij = exp(-C_i * dist);
  // Eq. (5)
  Real dRs_ij = eta_i * (df_ij - f_ij * C_i) * expCrij;
  // Eq. (19)
  Real d2Rs_ij = eta_i * (d2f_ij - 2 * df_ij * C_i + f_ij * C_i * C_i) * expCrij;
  //**********************************
  // If atom 1 is a polar H+, accumulate the derivative dR.
  if (topo->atomTypes[type1].mySCPISM_T->isHbonded == PH &&
      topo->atomTypes[type2].mySCPISM_T->isHbonded == PA &&
      (topo->atoms[atom1].residue_seq !=
       topo->atoms[atom2].residue_seq)) { // Polar
    Real E_i = 0.80; // Currently all polar H+ have this value
    Real g_i;
    if (topo->atoms[atom1].name == "HN" &&
        topo->atoms[atom2].name == "O")
      g_i = -0.378;
    else
      g_i = topo->atomTypes[type1].mySCPISM_T->g_i;
    Real g_j = topo->atomTypes[type2].mySCPISM_T->g_i;
    //topo->atoms[atom1].mySCPISM->polarFrac += g_i * g_j * f_ij * exp(
    //  -E_i * dist);
    Real expErij = exp(-E_i * dist);
    //**********************************
    // Eq. (7)
    dRs_ij +=  g_i * g_j * (df_ij - f_ij * E_i) * expErij;
    // Eq. (21)
    d2Rs_ij += g_i * g_j * (d2f_ij - 2 * df_ij * E_i + f_ij * E_i * E_i) * expErij;
    //**********************************
  }
  //Atom 2 Born Radius first derivative
  //**********************************
  Real C_j = topo->atomTypes[type2].mySCPISM_T->C_i;
  Real eta_j = topo->atoms[atom2].mySCPISM_A->eta;
  Real expCrji = exp(-C_j * dist);
  // Eq. (5)
  Real dRs_ji = eta_j * (df_ij - f_ij * C_j) * expCrji;
  // Eq. (19)
  Real d2Rs_ji = eta_j * (d2f_ij - 2 * df_ij * C_j + f_ij * C_j * C_j) * expCrji;
  //**********************************
  // If atom 2 is a polar H+, accumulate polar
  // fraction and derivative
  if (topo->atomTypes[type2].mySCPISM_T->isHbonded == PH &&
      topo->atomTypes[type1].mySCPISM_T->isHbonded == PA &&
      (topo->atoms[atom2].residue_seq !=
       topo->atoms[atom1].residue_seq)) {
    Real E_i = 0.80; // Currently all polar H+ have this value
    Real g_i;
    if (topo->atoms[atom2].name == "HN" && topo->atoms[atom1].name == "O")
      g_i = -0.378;
    else g_i = topo->atomTypes[type2].mySCPISM_T->g_i;
    Real g_j = topo->atomTypes[type1].mySCPISM_T->g_i;
    Real expErji = exp(-E_i * dist);
    //**********************************
    // Eq. (7)
    dRs_ji += g_i * g_j * (df_ij - f_ij * E_i) * expErji;
    // Eq. (21)
    d2Rs_ji += g_i * g_j * (d2f_ij - 2 * df_ij * E_i + f_ij * E_i * E_i) * expErji;
    //**********************************
  }      
  //Born radius
  Real Rs_i = topo->atoms[atom1].mySCPISM_A->bornRadius;
  Real rRs_i = 1.0 / Rs_i;
  Real Rs_j = topo->atoms[atom2].mySCPISM_A->bornRadius;
  Real rRs_j = 1.0 / Rs_j;
  //Screening
  Real D = dielecConst;
  Real K = (D - 1.0) / 2;
  Real alpha_i = topo->atomTypes[type1].mySCPISM_T->alpha;
  Real Ds_i = topo->atoms[atom1].mySCPISM_A->D_s;
  if(Ds_i == 0.0){
    // Eq. (9)
    Ds_i = (1.0 + D) / (1 + K * exp(-alpha_i * Rs_i)) - 1.0;
    topo->atoms[atom1].mySCPISM_A->D_s = Ds_i;
  }
  //Real Ds_i = topo->atoms[atom1].mySCPISM->D_s;
  Real rDs_i = 1.0 / Ds_i;
  // Eq. (12)
  Real dDs_i = alpha_i * (1.0 / (1.0 + D)) * (1 + Ds_i) * (D - Ds_i);
  // Eq. (18)
  Real d2Ds_i = alpha_i * (1.0 / (1.0 + D)) * (D - 1 - 2*Ds_i) * dDs_i;
  //Real expaRs = exp(-alpha_i * Rs_i);
  //Real dDs_i = ((1.0 + D) * K * alpha_i * expaRs) / ( (1 + K * expaRs) * (1 + K * expaRs) );
  //Real d2Ds_i = (2 * (1.0 + D) * K * K * alpha_i * alpha_i * expaRs * expaRs) / ( (1 + K * expaRs) * (1 + K * expaRs) * (1 + K * expaRs) )
  //               - ((1.0 + D) * K * alpha_i * alpha_i * expaRs) / ( (1 + K * expaRs) * (1 + K * expaRs) );

  Real alpha_j = topo->atomTypes[type2].mySCPISM_T->alpha;
  Real Ds_j = topo->atoms[atom2].mySCPISM_A->D_s;
  if(Ds_j == 0.0){    
    // Eq. (9)
    Ds_j = (1.0 + D) / (1 + K * exp(-alpha_j * Rs_j)) - 1.0;
    topo->atoms[atom2].mySCPISM_A->D_s = Ds_j;
  }
  //Real Ds_j = topo->atoms[atom2].mySCPISM->D_s;
  Real rDs_j = 1.0 / Ds_j;  
  // Eq. (12)
  Real dDs_j = alpha_j * (1.0 / (1.0 + D)) * (1 + Ds_j) * (D - Ds_j);  
  // Eq. (18)
  Real d2Ds_j = alpha_j * (1.0 / (1.0 + D)) * (D - 1 - 2*Ds_j) * dDs_j;
  //Real expaRs_j = exp(-alpha_j * Rs_j);
  //Real dDs_j = ((1.0 + D) * K * alpha_j* expaRs_j) / ( (1 + K * expaRs_j) * (1 + K * expaRs_j) );
  //Real d2Ds_j = (2 * (1.0 + D) * K * K * alpha_j * alpha_j * expaRs_j * expaRs_j) / ( (1 + K * expaRs_j) * (1 + K * expaRs_j) * (1 + K * expaRs_j) )
   //              - ((1.0 + D) * K * alpha_j * alpha_j * expaRs_j) / ( (1 + K * expaRs_j) * (1 + K * expaRs_j) );
  //Find derivatives of Energy
  Real q_i = topo->atoms[atom1].scaledCharge;
  Real q_j = topo->atoms[atom2].scaledCharge;
  // Eq. (10)
  Real dEoverR_ij = ( (0.5 * q_i * q_i * rRs_i * rRs_i *
                      (1.0 - rDs_i - Rs_i * rDs_i * rDs_i * dDs_i) * dRs_ij)
                      + (0.5 * q_j * q_j * rRs_j * rRs_j *
                          (1.0 - rDs_j - Rs_j * rDs_j * rDs_j * dDs_j) * dRs_ji) ) * rDist;
  // Eq. (16)
  Real d2E = (q_i * q_i * rRs_i) * ( 
                (dRs_ij * dRs_ij) * (dDs_i * dDs_i * rDs_i * rDs_i * rDs_i 
                                          - 0.5 * d2Ds_i * rDs_i * rDs_i 
                                             + dDs_i * rDs_i * rDs_i * rRs_i + rDs_i * rRs_i * rRs_i - rRs_i * rRs_i)
                  -0.5 * d2Rs_ij * (dDs_i * rDs_i * rDs_i + rDs_i * rRs_i - rRs_i)
               ) +
             (q_j * q_j * rRs_j) * ( 
                (dRs_ji * dRs_ji) * (dDs_j * dDs_j * rDs_j * rDs_j * rDs_j 
                                          - 0.5 * d2Ds_j * rDs_j * rDs_j 
                                             + dDs_j * rDs_j * rDs_j * rRs_j + rDs_j * rRs_j * rRs_j - rRs_j * rRs_j)
                  -0.5 * d2Rs_ji * (dDs_j * rDs_j * rDs_j + rDs_j * rRs_j - rRs_j)
               );

  //generate reduced hessian
  Matrix3By3 vec_rij_ij(rij, rij);   // outer products of vectors r_ij
  vec_rij_ij *= rDist * rDist;  //normalize

  Matrix3By3 I(1, 0, 0, 0, 1, 0, 0, 0, 1);   // now I is identity matrix
  // Eq. (22)
  Matrix3By3 H(I * dEoverR_ij + vec_rij_ij * (d2E - dEoverR_ij));

  return H;
}

