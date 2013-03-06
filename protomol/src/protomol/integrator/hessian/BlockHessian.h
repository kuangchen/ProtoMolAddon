/* -*- c++ -*- */
#ifndef BLOCKHESSIAN_H
#define BLOCKHESSIAN_H

#include <protomol/integrator/hessian/Hessian.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/BlockMatrix.h>

namespace ProtoMol {
  /**
   *
   * Calculates the Resudue Hessians or mass re-weighted Resudue Hessians
   * for the current force field.
   *
   */
  class BlockHessian : public Hessian {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  enum {MAX_ATOMS_PER_RES = 30};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    BlockHessian();
    ~BlockHessian();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class BlockHessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void clearBlocks();
    void initialResidueData(const GenericTopology *myTopo, int res_per_block, bool fullElect);
    void evaluateResidues(const Vector3DBlock *myPositions,
                            GenericTopology *myTopo,
                              bool simuMin);
    void evaluateNumericalResidues(const Vector3DBlock *myPositions,
                          const GenericTopology *myTopo);
    void evaluateBlockForces(const unsigned int blockStart, const unsigned int blockMax,
                                const Vector3DBlock *myPositions,
                                    const GenericTopology *myTopo, Vector3DBlock *blockForces);
    void evaluateBlocks(const Real cutoffDistance, 
                        const Vector3DBlock *myPositions,
                        GenericTopology *myTopo);
    void outputBlocks(const unsigned int ii, 
                           const unsigned int kk, const Matrix3By3 rhd);
    void evaluateInterBlocks(const Vector3DBlock *myPositions,
                             GenericTopology *myTopo);
    void outputMatrix(int i, int j, Real sqrtMassi, Real sqrtMassj, 
                        Matrix3By3 rha, int arrSz, double *basePoint);
  private:
    Real dihedralAngle(const int *aout, const Vector3DBlock &myPositions) const;
    void outputTorsions(const std::vector<Torsion> &torsions, const Vector3DBlock &myPositions);  

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends of class BlockHessian
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    unsigned int rsz; //size
    Real *sqrtMass;

  public:
    //Residues
    int num_residues, hess_array_size, hess_eig_size, *hess_array_point, *hess_eig_point;
    int *residues, *residues_max, *residues_alpha_c, 
        *residues_phi_n, *residues_psi_c;
    int *atom_residue, *atom_res_num;
    //Blocks
    int num_blocks, rpb, *blocks_max, *atom_block, *atom_block_num;
    //
    //Block Hessian
    vector<BlockMatrix> blocks;
    vector<BlockMatrix> adj_blocks;
    vector<BlockMatrix> non_adj_bond_blocks;
    vector<int> non_adj_bond_index;
    vector<BlockMatrix> adj_nonbond_blocks;
    vector<int> adj_nonbond_index;
    BlockMatrix electroStatics;
    bool fullElectrostatics;
    int memory_base, memory_blocks;

  };
}
#endif
