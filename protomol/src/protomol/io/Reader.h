/*  -*- c++ -*-  */
#ifndef READER_H
#define READER_H

#include <protomol/io/File.h>

namespace ProtoMol {
  //_________________________________________________________________ Reader
  /**
   * Base class of readers
   */
  class Reader : public File {
  protected:
    Reader();
    Reader(const std::string &filename);
    /// To open with special file flags, std::ios::in is set
    Reader(std::ios::openmode mode);
    /// To open with special file flags, std::ios::in is set
    Reader(std::ios::openmode mode, const std::string &filename);

  public:
    /// Simple test, true if the format might be correct/readable
    virtual bool tryFormat() = 0;
    virtual bool read() = 0;
  };
}

#endif /* READER_H */
