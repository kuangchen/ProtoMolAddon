/*  -*- c++ -*-  */
#ifndef FILE_H
#define FILE_H

#include <string>
#include <vector>
#include <protomol/base/Exception.h>

#ifdef HAVE_LIBFAH
#include <cbang/os/File.h>
#else
#include <fstream>
#endif

namespace ProtoMol {
  //_________________________________________________________________ File
  /**
   * Abstract base class for all I/O; readers and writers. The readers and
   * writes are intend to act STL alike to stream into or from a supported
   * structure or container.
   *
   * NB:
   * - New writer or reader should never inherit directly from File, but
   *   from Reader or Writer.
   * - Reading binaries one should always use read(), rather directly
   *   file.read, since some compilers like Sun WorkShop CC have problems.
   * - File objects can be used as ios_base objects inside expression
   *   (e.g., while(dcdReader >> xyz){ ... } )
   */
  class File {
  protected:
    std::ios::openmode mode;
    std::string filename;
    std::string comment;
#ifdef HAVE_LIBFAH
    cb::File file;
#else
    std::fstream file;
#endif

    File();
    File(std::ios::openmode mode);
    File(std::ios::openmode mode, const std::string &filename);

  public:
    virtual ~File();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class File
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void setFilename(const std::string &filename) {this->filename = filename;}
    std::string getFilename() const {return filename;}
    const std::string &getComment() const {return comment;}
    void setComment(const std::string &comment) {this->comment = comment;}

    bool isAccessible();
    bool open(const std::string &filename);
    bool open(const std::string &filename, std::ios::openmode mode);
    bool open();
    void close();
    bool is_open();
    bool eof() {return file.eof();}

    // enable expression testing
    operator void*() const {return !*this ? 0 : const_cast<File *>(this);}
    bool operator!() const {return file.fail();}

    void write(const char *c, std::streamsize count);
    // Redirect of fstream::read  (Sun WorkShop CC does not properly read
    // more than one char ...)
    void read(char *c, std::streamsize count);
    std::string getline();
    unsigned int getLineTokens(std::vector<std::string> &tokens);
  };
}

#endif /* FILE_H */
