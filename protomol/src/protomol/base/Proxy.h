/*  -*- c++ -*-  */
#ifndef PROXY_H
#define PROXY_H

namespace ProtoMol {
//_________________________________________________________________ Proxy
/**
 * Book keeping for distribution and reduction of data structures in
 * parallel environment.
 */
  class Proxy {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Proxy() : myCount(0) {}
    virtual ~Proxy() {}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Proxy
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void distribute()  const {myCount++;}
    void reduce()      const {myCount--;}
    bool distributed() const {return myCount > 0;}
    int distCount()       const {return myCount;}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    mutable int myCount;
  };
}
#endif /* PROXY_H */
