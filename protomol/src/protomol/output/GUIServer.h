#ifndef PROTOMOL_GUI_SERVER_H
#define PROTOMOL_GUI_SERVER_H

#define COMM_MAGIC   0xF01DBAAD
#define COMM_VERSION 2
#define COMM_PORT    0xCE11
#define COMM_NNCOORD 050363
#define NUMCONN      20

#include "stdtypes.h"
#include <sys/types.h>

#if defined _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock.h>
typedef int socklen_t;  //  Unix socket length

#else
typedef int SOCKET;
#include <pthread.h>
#endif

namespace ProtoMol {

  struct FAH_ATOM {
    uint8_t type[4];
    float charge;
    float radius;
#if defined _WIN32
  };
#else
  } __attribute__((packed));
#endif

  struct FAH_BOND {
    uint32_t a; /* rule: a < b */
    uint32_t b;
#if defined _WIN32
  };
#else
  } __attribute__((packed));
#endif

  struct FAH_XYZ{
    float x;
    float y;
    float z;
#if defined _WIN32
  };
#else
  } __attribute__((packed));
#endif

  struct FAH_INFO{
    uint32_t magic;
    uint32_t version;
    uint8_t  name[64];
    int64_t  timestamp;
    uint64_t iterations;
    uint32_t frames;
    uint32_t atom_count;
    uint32_t bond_count;
#if defined _WIN32
  };
#else
  } __attribute__((packed));
#endif

  struct FAH_CURRENT{
    uint32_t magic;
    uint32_t version;
    int64_t  timestamp;
    uint64_t iterations_done;
    uint32_t frames_done;
    float    energy;
    float    temperature;
#if defined _WIN32
  };
#else
  } __attribute__((packed));
#endif

  class GUIServer {
  public:
    typedef enum {
      GS_NO_REQUEST,
      GS_META_REQUEST,
      GS_COORD_REQUEST,
    } request_t;

  protected:
#ifdef _WIN32
    HANDLE thread;
    HANDLE mutex;
#else
    pthread_t thread;
    pthread_mutex_t mutex;
#endif

    request_t request;
    int shutdown;
    SOCKET connectlist[NUMCONN];
    bool connectlist_sent_curr[NUMCONN];
    fd_set socks;
    unsigned short com_port;
    unsigned int com_port_range;
    uint32_t no_new_coord;

  public:
    FAH_INFO info;
    FAH_CURRENT current;
    FAH_ATOM *atoms;
    FAH_BOND *bonds;
    FAH_XYZ *xyz;

    GUIServer(const char *name, int natoms, int nbonds, int port, int prange);
    ~GUIServer();

    void startUpdate();
    void endUpdate();
    request_t getRequest() {return request;}
    void startServer();

  protected:
    int sends(SOCKET socket, char *data, int length);
    void lock();
    void unlock();
    int sendMetadata(SOCKET socket);
    int sendCoordinates(SOCKET socket);
    int sendNoNewCoord(SOCKET socket);
    void serverThread();
#ifdef _WIN32
    static DWORD WINAPI callServerThread(LPVOID);
#else
    static void *callServerThread(void *param);
#endif
  };
}

#endif //  PROTOMOL_GUI_SERVER_H
