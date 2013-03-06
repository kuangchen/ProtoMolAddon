/*  -*- c++ -*-  */
#ifndef CELLLISTENUMERATOR_H
#define CELLLISTENUMERATOR_H

namespace ProtoMol {
  /**
   * The cell enumerator implements an iterator over the cell list
   * in respect to the give cutoff. It iterates over all
   * possible cell pairs such that all interactions are
   * considered of at least the cutoff.
   */

  struct CellPair {int first; int second;};

  template<class TBoundaryConditions, class TCellManager>
  class CellListEnumerator {};

}
#endif /* CELLLISTENUMERATOR_H */
