/*  -*- c++ -*-  */
#ifndef PROTOMOL_OUTPUT_COLLECTION_H
#define PROTOMOL_OUTPUT_COLLECTION_H

#include <list>

namespace ProtoMol {
  class Output;
  class OutputFactory;
  class ProtoMolApp;

  // / Container class for Output objects invoked at application level.
  class OutputCollection  {
    friend class OutputFactory;

    typedef std::list<Output *> Container;
    Container outputList;

    const ProtoMolApp *app;

  public:
    OutputCollection() : app(0) {}
    ~OutputCollection();

    // / Initialize all Output object
    void initialize(const ProtoMolApp *app);

    // / Invoke all Output objects with run().  Returns true if an Output ran.
    bool run(int step);

    // / Finalize all Outout object
    void finalize(int step);

    // / Add new Output object to the collection
    int getNext() const;

    void adoptOutput(Output *output);

    // / Iterators, const
    typedef Container::const_iterator const_iterator;
    const_iterator begin() const {return outputList.begin();}
    const_iterator end()   const {return outputList.end();}

  private:
    // / Iterator
    typedef Container::iterator iterator;
    iterator begin()       {return outputList.begin();}
    iterator end()         {return outputList.end();}
  };
}
#endif //  PROTOMOL_OUTPUT_COLLECTION_H
