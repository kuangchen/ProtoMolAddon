#include <protomol/base/SystemUtilities.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h> // for GetFullPathName

//____ Define the missing symbols from <unistd.h> for M$ ....
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#define access _access
#include <fcntl.h>
#include <io.h>
#define F_OK 0
#define W_OK 2
#define R_OK 4

#define S_ISDIR(x) (((x) & _S_IFDIR) != 0)

#else // _WIN32

#include <libgen.h>
#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif

#include <fstream>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_LIBFAH
#include <fah/checksum/overrides.h>
#endif

using namespace std;

namespace ProtoMol {
  namespace SystemUtilities {
#ifdef _WIN32
    const string path_separators  = "/\\";
    const char path_separator  = '\\';
    const char path_delimiter  = ';';

#else
    const string path_separators = "/";
    const char path_separator = '/';
    const char path_delimiter = ':';

#endif

    string basename(const string &path) {
      string::size_type pos = path.find_last_of(path_separators);

      if (pos == string::npos) return path;
      else return path.substr(pos + 1);
    }


    string dirname(const string &path) {
      string::size_type pos = path.find_last_of(path_separators);

      if (pos == path.length() - 1)
        pos = path.find_last_of(path_separators, pos - 1);

      if (pos == string::npos) return ".";
      else if (pos == 0) return "/";
      else return path.substr(0, pos);
    }


    bool exists(const string &path) {
      struct stat buf;
      return stat(path.c_str(), &buf) != -1;
    }


    void mkdir(const string &path, bool withParents) {
      if (withParents) {
        string parent = dirname(path);

        if (parent != ".") {
          if (!isDirectory(parent)) {
            if (exists(parent))
              THROWS("'" << parent << "' exists but is not a directory");

            mkdir(parent, true);
          }
        }
      }

#ifdef _WIN32
      if (::mkdir(path.c_str()) == -1)
#else
      if (::mkdir(path.c_str(), 0777) == -1)
#endif
        THROWS("Failed to create directory '" << path << "'");
    }


    bool isDirectory(const string &path) {
      struct stat buf;
      if (stat(path.c_str(), &buf) == -1) return false;
      return S_ISDIR(buf.st_mode);
    }


    bool ensureDirectory(const string &path) {
      if (!isDirectory(path)) {
        if (exists(path))
          THROWS("'" << path << "' exists but is not a directory");

        mkdir(path, true);

        return true;
      }

      return false;
    }

    bool chdir(const string &path) {
#ifdef _WIN32
      if (::_chdir(path.c_str()) < 0)
#else
      if (::chdir(path.c_str()) < 0)
#endif
        THROWS("chdir(" << path << ") failed");

      return true;
    }


    bool isAccessible(const string &fileName) {
      return ::access(fileName.c_str(), F_OK) == 0;
    }


    void splitFileName(const string &filename, string &dirname,
                       string &basename, string &extension) {
      string::size_type pos = filename.rfind(PATHSEP);

      if (pos == string::npos) basename = filename;
      else {
        dirname = filename.substr(0, pos);
        basename = filename.substr(pos + 1);
      }

      pos = basename.rfind('.');

      if (pos != string::npos) {
        extension = basename.substr(pos + 1);
        basename = basename.substr(0, pos);
      }
    }


    unsigned int getFileSize(const string &filename) {
      ifstream f(filename.c_str());
      if (!f.is_open()) THROWS("Failed to open '" << filename << "'");

      f.seekg(0, ios::end);
      return (unsigned int)f.tellg();
    }
  }


  static void (*myAbortFunction)() = NULL;

  void protomolAbort() {
    if (myAbortFunction) (*myAbortFunction)();

    THROW("ABORT");
  }


  void setProtomolAbort(void (*abortFunction)()) {
    myAbortFunction = abortFunction;
  }


  static void (*myExitFunction)() = NULL;

  void protomolExit() {
    if (myExitFunction) (*myExitFunction)();

    THROW("EXIT");
  }


  void setProtomolExit(void (*exitFunction)()) {
    myExitFunction = exitFunction;
  }


  static void (*myStartSerial)(bool) = NULL;

  void protomolStartSerial(bool exludeMaster) {
    if (myStartSerial != NULL)
      (*myStartSerial)(exludeMaster);
  }


  void setProtomolStartSerial(void (*startSerialFunction)(bool)) {
    myStartSerial = startSerialFunction;
  }


  static void (*myEndSerial)(bool) = NULL;


  void protomolEndSerial(bool exludeMaster) {
    if (myEndSerial != NULL)
      (*myEndSerial)(exludeMaster);
  }


  void setProtomolEndSerial(void (*endSerialFunction)(bool)) {
    myEndSerial = endSerialFunction;
  }


  struct Endian {
    // Helper class to make sure that we get endianess correct ... M$
    static bool isLittleEndian() {
      unsigned int tmp = 1;
      return 0 != *(reinterpret_cast<const char *> ( &tmp));
    }
  };
  const bool ISLITTLEENDIAN = Endian::isLittleEndian();


  string getCanonicalPath(const string &path) {
    char buf[4096];

#ifdef _WIN32
    char *finalpart;
    DWORD len = GetFullPathName(path.c_str(), 4096, buf, &finalpart);
    if (len == 0 || len > 4095)
      THROW(string("GetFullPathName '") + path + "' failed.");

    return buf;

#else
    char tmp[path.length() + 3];

    // The file might not exist yet but its directory must.
    strcpy(tmp, path.c_str());
    string dir = dirname(tmp);

    if (!realpath(dir.c_str(), buf))
      THROW(string("realpath '") + path + "' failed.");

    strcpy(tmp, path.c_str());

    return string(buf) + "/" + basename(tmp);
#endif
  }


  bool SystemUtilities::unlink(const string &path) {
#ifdef HAVE_LIBFAH
    return fah_unlink(path.c_str()) == 0;
#else
    return ::unlink(path.c_str()) == 0;
#endif
  }


  void SystemUtilities::rename(const string &src, const string &dst) {
    unlink(dst);
#ifdef HAVE_LIBFAH
    int error = fah_rename(src.c_str(), dst.c_str());
#else
    int error = ::rename(src.c_str(), dst.c_str());
#endif

    if (error)
      THROWS("Failed to rename '" << src << "' to '" << dst << "': "
             << strerror(errno));
  }
}
