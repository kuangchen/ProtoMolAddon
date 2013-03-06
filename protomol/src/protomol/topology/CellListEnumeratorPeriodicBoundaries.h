/* -*- c++ -*- */
#ifndef CELLLISTENUMERATOR_PERDIOCBOUNDARIES_H
#define CELLLISTENUMERATOR_PERDIOCBOUNDARIES_H

#include <protomol/topology/CellListEnumerator.h>

#include <protomol/topology/CubicCellManager.h>
#include <protomol/topology/PeriodicBoundaryConditions.h>
#include <protomol/topology/Topology.h>
#include <set>
namespace ProtoMol {
  //_______________________________________ CellListEnumeratorPeriodicBoundaries

  /**
   * Specialization of the cell enumerator for periodic boundary conditions
   * and cubic cell manager
   */
  template<>
  class CellListEnumerator<PeriodicBoundaryConditions, CubicCellManager> {
  public:
    CellListEnumerator() : myCutoff(-1.0),
      myCellSize(Vector3D(-1.0, -1.0,
                          -1.0)), myMax(CubicCellManager::Cell(-1, -1, -1)) {}

    void initialize(const Topology<PeriodicBoundaryConditions,
                                   CubicCellManager> *topo,
                    Real cutoff) {
      the_beginning = topo->cellLists.begin();
      the_end = topo->cellLists.end();

      i = the_beginning;
      j = i;

      the_first = 0;
      the_counter = the_first;
      myCellListStruct = &(topo->cellLists);

      if (myCutoff != cutoff ||
          topo->cellManager.getCellSizeVector() != myCellSize ||
          !(myMax == topo->cellManager.findCell(topo->max - topo->min))) {
        myCellSize = topo->cellManager.getCellSizeVector();
        myCutoff = cutoff;

        myMax = topo->cellManager.findCell(topo->max - topo->min);

        Real cutoff2 = cutoff * cutoff;

        CubicCellManager::Cell zero(0, 0, 0);
        myDeltaList.clear();
        int nx = (int)(cutoff / myCellSize.c[0] + 1.0 + Constant::EPSILON);
        int ny = (int)(cutoff / myCellSize.c[1] + 1.0 + Constant::EPSILON);
        int nz = (int)(cutoff / myCellSize.c[2] + 1.0 + Constant::EPSILON);
        std::set<CubicCellManager::Cell> tmpDelta;
        tmpDelta.insert(zero);
        // Do not consider deltas bigger than the dimesion of the simulation box
        int n0 = std::min(nx, topo->cellLists.getDimX());
        int n1 = std::min(ny, topo->cellLists.getDimY());
        int n2 = std::min(nz, topo->cellLists.getDimZ());
        Real xx = myCellSize.c[0] * myCellSize.c[0];
        Real yy = myCellSize.c[1] * myCellSize.c[1];
        Real zz = myCellSize.c[2] * myCellSize.c[2];
        for (int k = -n0; k <= n0; k++) {
          int x = abs(k) - 1; if (x < 0) x = 0;
          Real d0 = x * x * xx;
          for (int l = -n1; l <= n1; l++) {
            int y = abs(l) - 1; if (y < 0) y = 0;
            Real d1 = d0 + y * y * yy;
            for (int m = -n2; m <= n2; m++) {
              int z = abs(m) - 1; if (z < 0) z = 0;
              if (d1 + z * z * zz < cutoff2) {
                CubicCellManager::Cell delta(k, l, m);
                CubicCellManager::Cell minimalDelta =
                  myCellListStruct->basisCell(delta);
                if (tmpDelta.find(minimalDelta) == tmpDelta.end()) {
                  myDeltaList.push_back(minimalDelta);
                  tmpDelta.insert(minimalDelta);
                }
              }
            }
          }
        }

        std::sort(myDeltaList.begin(), myDeltaList.end());
      }
      the_last = myDeltaList.size();
    }

    /// retrieve the current cell pair
    void get(CellPair &cp) {cp.first = i->second; cp.second = j->second;}
    /// retrieve the current cell pair
    void get(int &a, int &b) {a = i->second; b = j->second;}
    /// if the cells of the current cell pair are the same
    bool notSameCell() {return i != j;}
    /// reached the end of the list of cell pairs
    bool done() {return i == the_end;}
    /// goto the end of pair of cells with the same first cell
    void gotoEndPair() {j = the_end;};

    /// advance by inc in the cell list in respect to first cell
    void nextNewPair(int inc) {
      if (inc < 1)
        return;
      while (inc > 0 && i != the_end) {
        ++i;
        --inc;
        while (i != the_end && i->second < 0)
          ++i;
      }

      j = i;
    }

    /// get next pair
    void next(void) {
      do
        if (i != the_end) {
          j = the_end;
          while (the_counter != the_last &&
                 the_end ==
                 (j =
                    myCellListStruct->findPeriodic(i->first +
                      myDeltaList[the_counter++]))) ;

          if (j == the_end) {
            ++i;
            while (i != the_end && i->second < 0)
              ++i;

            j = i;
            the_counter = the_first;
          }
        }
      while (i > j);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    CubicCellManager::CellListStructure::const_iterator i, j, the_beginning,
                                                        the_end;
    Real myCutoff;
    std::vector<CubicCellManager::Cell> myDeltaList;
    int the_first, the_last, the_counter;
    const CubicCellManager::CellListStructure *myCellListStruct;
    Vector3D myCellSize;
    CubicCellManager::Cell myMax;
  };
}
#endif /* CELLLISTENUMERATOR_PERDIOCBOUNDARIES_H */
