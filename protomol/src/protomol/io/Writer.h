/*  -*- c++ -*-  */
#ifndef WRITER_H
#define WRITER_H

#include <protomol/io/File.h>

namespace ProtoMol {
  //_________________________________________________________________ Writer
  /**
   * Base clas of writers
   */
  class Writer : public File {
  protected:
    Writer();
    Writer(const std::string &filename);
    /// To open with special file flags, std::ios::out is set
    Writer(std::ios::openmode mode);
    /// To open with special file flags, std::ios::out is set
    Writer(std::ios::openmode mode, const std::string &filename);
  };
}

#endif /* WRITER_H */
