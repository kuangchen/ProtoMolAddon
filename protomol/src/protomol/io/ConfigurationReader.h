/*  -*- c++ -*-  */
#ifndef CONFIGURATIONREADER_H
#define CONFIGURATIONREADER_H

#include <protomol/io/Reader.h>
#include <protomol/config/Configuration.h>

namespace ProtoMol {
  //_________________________________________________________ConfigurationReader
  /*
   * Reads and parses a ProtoMol configuaration file. The acutal parsing is
   * delegated to each entry of the configuration container. The entries have
   * a keyowrd (identifier) and a value with associated type and constraint.
   * The parsing is implemented in the traits of the types.
   */
  class ConfigurationReader : public Reader {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors (both default here), assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ConfigurationReader();
    ConfigurationReader(const std::string &filename);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool open() {return File::open();}
    virtual bool open(const std::string &filename)
    {return File::open(filename);}
    virtual bool open(const char *filename) {return File::open(filename);}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Reader
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool tryFormat();
    virtual bool read() {return false;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Configuration
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    bool read(Configuration &config);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Friends
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    friend ConfigurationReader &operator>>(ConfigurationReader &configReader,
                                           Configuration &config);
  };
}
#endif /* CONFIGURATIONXYZREADER_H */
