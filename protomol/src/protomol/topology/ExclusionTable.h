/* -*- c++ -*- */
#ifndef EXCLUSIONTABLE_H
#define EXCLUSIONTABLE_H

#include <vector>
#include <map>
#include <set>

#include <protomol/base/PMConstants.h>
#include <protomol/type/SimpleTypes.h>

namespace ProtoMol {
  //________________________________________ ExclusionClass
  /**
   * Defines exclusions between two atoms, intra-molecular
   * They are used to implement (bulding blocks) the exclusion type of
   * system: ONE2, ONE3, ONE4 and ONE4MODIFIED
   */
  enum ExclusionClass {
    EXCLUSION_NONE = 0,     ///< do not exclude
    EXCLUSION_MODIFIED = 1, ///< do not exlcude, but modify
    EXCLUSION_FULL = 2      ///< do exclude
  };

  //________________________________________ ExclusionPair
  /**
   * Struct to keep track and sort all kind exclusions
   */
  struct ExclusionPair {
    ExclusionPair();
    ExclusionPair(int a, int b);
    ExclusionPair(int a, int b, ExclusionClass c);
    bool operator<(const ExclusionPair &e) const;

    int a1, a2;
    ExclusionClass excl;
  };

  //________________________________________ ExclusionTable
  /**
   * Defines the table of exclusions. ExclusionTable assumes that
   * most exclusions have a difference (delta) of the indices's of the two atoms
   * less than Constant::FASTDELTAMAX, it keeps also track of the maximum
   * difference. If the distance is less than Constant::FASTDELTAMAX a simple
   * array look-up table is used, otherwise a map will be used.
   */
  class ExclusionTable {
  public:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ExclusionTable();
    ~ExclusionTable();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ExclusionTable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /**
     * Set the number of atoms in the system (not the number of exclusions).
     * This invalidates the data, so a clear must be performed after this 
     * before the table is used again.
     */
    void resize(int count);

    /// Clear all exclusions from the table.
    void clear();


    /// Add a new exclusion between atom1 and atom2 of the given type.
    void add(int atom1, int atom2, ExclusionClass type);

    /// Check for an exclusion between atom1 and atom2, and return its type.
    ExclusionClass check(int atom1, int atom2) const;

    /**
     * Check if there could be an exclusion between atom1 and atom2.  If
     * this function returns false, there is no exclusion, but a true
     * result doesn't give any information.  Requires that atom1 $<$ atom2.
     */
    bool checkReallyFast(int atom1, int atom2) const;

    /// Returns true if there are no exclusions.
    bool empty() const {return myTable.empty();}

    const std::vector<ExclusionPair> &getTable() const {return myTable;}

    /// Clean up temporaries and optimize space, to be called after the
    /// exclusion table is build
    void optimize();

    /// Maximal distance between two excluded atom pair
    int getMaxDelta() const {return myMaxDelta;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    ExclusionClass *lowDeltas;
    char *anyExclsForThisDelta;
    std::map<PairInt, ExclusionClass> highDeltas;
    int myMaxDelta;
    std::vector<ExclusionPair> myTable;
    std::set<ExclusionPair> mySet;
    int myCurrentSize;
    int myFastDeltaMax;
  };
  //________________________________________ INLINES

  inline bool ExclusionTable::checkReallyFast(int atom1, int atom2) const {
    int delta = atom2 - atom1;
    if (delta > myMaxDelta) return false;
    return true;
  }

  inline ExclusionClass ExclusionTable::check(int atom1, int atom2) const {
    if (atom2 < atom1)
      std::swap(atom1, atom2);
    int delta = atom2 - atom1;
    if (delta > myMaxDelta) return EXCLUSION_NONE;
    if (delta < myFastDeltaMax)
      return lowDeltas[atom1 * myFastDeltaMax + delta];
    else {
      if (anyExclsForThisDelta[delta] < 1) return EXCLUSION_NONE;
      std::map<PairInt, ExclusionClass>::const_iterator i = highDeltas.find(
        PairInt(atom1, delta));
      if (i == highDeltas.end())
        return EXCLUSION_NONE;
      return i->second;
    }
  }
}
#endif /* EXCLUSION_H */
