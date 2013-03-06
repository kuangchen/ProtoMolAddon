#include "GUIServer.h"

#include <protomol/base/Exception.h>

#include <iostream>

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#ifndef _WIN32
#include <sched.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
// #include <errno.h>
#include <fcntl.h>
#endif

//  Redefine some types and constants based on OS
#ifndef _WIN32
#define INVALID_SOCKET -1        //  WinSock invalid socket
#define SOCKET_ERROR   -1        //  Basic WinSock error
#define closesocket(s) close(s); //  Unix uses file descriptors, WinSock doesn't
#endif

/* This code sets up a server. When the server gets a connection,
 * it reads a 4 byte request code. We then set a flag with the value of
 * the requested data. This flag is examined periodically by the
 * core and the core may choose to fill in values when it is good
 * and ready. When it does so, the flag gets reset indicating
 * that the data has been filled in. If the client can not wait,
 * then the core wrapper has no choice but to return old valu
 * since the science code may be busy. For now, I'm assuming the
 * client code can smoothly handle delays from the core. So we
 * just use one mutex that indicates a request has been made.
 * While this mutex is locked, the core wrapper code will not
 * touch the coordinates etc.
 */


using namespace std;
using namespace ProtoMol;// FAH;

GUIServer::GUIServer(const char *name, int natoms, int nbonds, int port,
                     int prange) : request(GS_COORD_REQUEST), shutdown(0) {

  /* Fill in stuff in that the science code need not bother with
   * Also, init to non-junk values, since the viewer may updat
   * before the science code has had a chance to fill in anything
   */
  info.version = COMM_VERSION;
  info.magic   = COMM_MAGIC;
  strncpy((char *)info.name, name, 64);
  info.name[63] = 0;
  info.timestamp = time(0);
  info.frames     = 0;
  info.iterations = 0;
  info.atom_count = natoms;
  info.bond_count = nbonds;

  current.energy      = 0.0f;
  current.temperature = 0.0f;
  current.magic   = COMM_MAGIC;
  current.version = COMM_VERSION;
  current.iterations_done = 0;
  current.frames_done = 0;
  current.timestamp = info.timestamp;

  atoms = new FAH_ATOM[natoms];
  bonds = new FAH_BOND[nbonds];
  xyz   = new FAH_XYZ[natoms];

  no_new_coord = COMM_NNCOORD;
  com_port = port;
  com_port_range = prange;
}


GUIServer::~GUIServer() {
  shutdown = 1;

  //  Wait for server thread
#ifdef _WIN32
  if (thread) WaitForSingleObject(thread, 5000); //  Upto 5
#endif

#ifdef _WIN32
  if (thread) CloseHandle(thread);
  thread = 0;
  if (mutex) CloseHandle(mutex);
  mutex = 0;

  WSACleanup();
#endif

  delete [] atoms;
  delete [] bonds;
  delete [] xyz;

  atoms = 0;
  bonds = 0;
  xyz = 0;
}


int GUIServer::sends(SOCKET socket, char *data, int length) {
  int sendRet;
  int numSent = 0;

  //  non-blocking so wait for buffer
  while (numSent < length){
    sendRet = ::send(socket, &data[numSent], length-numSent, 0);
    if (sendRet < 0) {
#ifndef WIN32
      if (errno != EWOULDBLOCK) return -1;
#else
      if (WSAGetLastError() != WSAEWOULDBLOCK) return -1;
#endif
    } else numSent += sendRet;
  }

  return 0;
}


int GUIServer::sendMetadata(SOCKET socket) {
  return
    sends(socket, (char *)&info, sizeof(FAH_INFO)) ||
    sends(socket, (char *)atoms, sizeof(FAH_ATOM) * info.atom_count) ||
    sends(socket, (char *)bonds, sizeof(FAH_BOND) * info.bond_count);
}


int GUIServer::sendCoordinates(SOCKET socket) {
  return
    sends(socket, (char *)&current, sizeof(FAH_CURRENT)) ||
    sends(socket, (char *)xyz, sizeof(FAH_XYZ) * info.atom_count);
}


int GUIServer::sendNoNewCoord(SOCKET socket) {
  return
    sends(socket, (char *)&no_new_coord, sizeof(uint32_t));
}


void GUIServer::serverThread() {
  struct sockaddr_in addr, r_addr;
  SOCKET lsocket = INVALID_SOCKET;
  SOCKET asocket = INVALID_SOCKET;
  socklen_t len = sizeof(r_addr);   //  The length of our remote addr

  if ((lsocket = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) == INVALID_SOCKET)
    THROW("Invalid socket!");

  //  Setup the local address variabl
  memset((void *)&addr, 0, sizeof(addr));
  addr.sin_family      = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port        = htons(com_port);// COMM_PORT);

  //  Name the local socket
  if (bind(lsocket, (struct sockaddr *)&addr, sizeof(addr)) == SOCKET_ERROR) {
    // try one port up if failed (use for dual/quad core machines)
    if (com_port_range > 1) {
      int bret = SOCKET_ERROR;
      for (unsigned int i = 1; i < com_port_range && bret == SOCKET_ERROR; i++){
        addr.sin_port = htons(com_port+i);
        bret = bind(lsocket, (struct sockaddr *)&addr, sizeof(addr));
      }

      if (bret == SOCKET_ERROR) {
        closesocket(lsocket);
        THROWS("Could not bind to local socket: " << strerror(errno));
      }

    } else {
      closesocket(lsocket);
      THROWS("Could not bind to local socket: " << strerror(errno));
    }
  }

  //  Set the socket to listen for a connection
  if (listen(lsocket, SOMAXCONN) == SOCKET_ERROR) {//
    closesocket(lsocket);
    THROWS("Could not listen to local socket: " << strerror(errno));
  }

  //  Wait for a connection
  // initialize connection list
  for (int listnum = 0; listnum < NUMCONN; listnum++) {
    connectlist[listnum] = INVALID_SOCKET;
    connectlist_sent_curr[listnum] = true;  // don't send until real dat
  }

  // ####Main server loo
  //  While not close connection
  while (!shutdown) {
    int highsock = lsocket;
    int close_connection = 0;
    //  Wait for a connection
    // initialize file descriptor
    FD_ZERO(&socks);
    FD_SET(lsocket,&socks);

    for (int listnum = 0; listnum < NUMCONN; listnum++){
      if (connectlist[listnum] != INVALID_SOCKET){
        FD_SET(connectlist[listnum],&socks);
        if (connectlist[listnum] > highsock)
          highsock = connectlist[listnum];
      }
    }

    // wait for socket activity
    int readsocks =
      select(highsock + 1, &socks, (fd_set *)0, (fd_set *)0, 0);
    if (readsocks < 0) THROWS("Socket non-blocking error: " << strerror(errno));
    if (readsocks) {
      // tivity on new connection?
      if (FD_ISSET(lsocket, &socks)) {
        int socket_error = 0;
        if ((asocket =
             accept(lsocket,
                    (struct sockaddr *)&r_addr, &len)) != INVALID_SOCKET){
          // t non blocking
#ifndef WIN32
          int opts;
          if ((opts = fcntl(asocket,F_GETFL)) >= 0) {
            opts = (opts | O_NONBLOCK);
            if (fcntl(asocket,F_SETFL,opts) >= 0) {
#else
              u_long on = 1;
              if (ioctlsocket(asocket, FIONBIO, &on) != SOCKET_ERROR) {
#endif
                // ut in list
                for (int listnum = 0; (listnum < NUMCONN) && (asocket != -1);
                     listnum ++) {
                  if (connectlist[listnum] == INVALID_SOCKET){
                    printf("Connection accepted:   FD=%d; Slot=%d\n", asocket,
                           listnum);
                    connectlist[listnum] = asocket;
                    asocket = -1;
                  }
                }

                if (asocket != -1) {
                  //  No room left in the queu
                  printf("No room left for new client.\n");
                  socket_error = 1;
                }
#ifndef WIN32
              } else {
                printf("Error on fcntl(F_SETFL)\n");
                socket_error = 1;
              }

            } else {
              printf("Error on fcntl(F_GETFL)\n");
              socket_error = 1;
            }
#else
          } else {
            printf("Error on ioctl\n");
            socket_error = 1;
          }
#endif
        } else {
          printf("Error on accept\n");
          socket_error = 1;
        }

        if (socket_error) closesocket(asocket);
      }

      //  handle data for all entries in queu
      unsigned recvrequest;
      for (int listnum = 0; listnum < NUMCONN; listnum++) {
        if (connectlist[listnum] != INVALID_SOCKET &&
            FD_ISSET(connectlist[listnum], &socks)) {
          close_connection = 0;

          if (recv(connectlist[listnum], (char *)&recvrequest, 4, 0) == 4) {
            //  Switch on request
            switch (recvrequest) {
            case 0: if (shutdown) break; //  send metadat
              if (sendMetadata(connectlist[listnum])) close_connection = 1;
              break;

            case 1: if (shutdown) break; //  send dat
              lock();
              //  Globally new daya? Required in case version 1 client in set
              if (request != GS_COORD_REQUEST)
                for (int clist = 0; clist < NUMCONN; clist++)
                  //  set data available to all client
                  connectlist_sent_curr[clist] = false;

              request = GS_COORD_REQUEST;

              if (sendCoordinates(connectlist[listnum])) close_connection = 3;

              unlock();
              break;

            case 2: if (shutdown) break; //  send data if new availabl
              lock();
              if (request == GS_COORD_REQUEST &&
                  connectlist_sent_curr[listnum]) {
                //  no new data? AND this connection has sent it
                unlock();
                if (sendNoNewCoord(connectlist[listnum])) close_connection = 2;

              } else {
                //  Globally new daya, or just this client?
                if (request != GS_COORD_REQUEST) //  global
                  for (int clist = 0; clist < NUMCONN; clist++)
                    //  set data available to all client
                    connectlist_sent_curr[clist] = false;

                //  service request
                request = GS_COORD_REQUEST;
                //  flag this client has had new dat
                connectlist_sent_curr[listnum] = true;

                if (sendCoordinates(connectlist[listnum])) close_connection = 3;
                unlock();
              }
              break;

            case 99: close_connection = 4; break;

            default: close_connection = 5;
            }

          } else close_connection = 6;

          if (close_connection) {
            printf("Connection closed:   FD=%d; Slot=%d; Err=%d;\n",
                   connectlist[listnum], listnum, close_connection);
            closesocket(connectlist[listnum]);
            connectlist[listnum] = INVALID_SOCKET;
          }
        }
      }
      //  end handle dat
    }
  }

  closesocket(lsocket);
  return;
}


#ifdef _WIN32
DWORD WINAPI GUIServer::callServerThread(LPVOID param)
#else
  void *GUIServer::callServerThread(void *param)
#endif
{
  GUIServer *server = dynamic_cast<GUIServer *>((GUIServer *)param);
  if (server){
    try {
      server->serverThread();
    } catch (const Exception &e) {
      cerr << e.getMessage();
    }
  }

  return 0;
}

void GUIServer::startServer() {
  shutdown = 0;

#ifdef _WIN32
  WSADATA wsa;
  DWORD tid;

  try {
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != NO_ERROR)
      THROW("Error starting winsock");

    if (LOBYTE(wsa.wVersion) != 2 || HIBYTE(wsa.wVersion) != 2)
      THROW("Error need winsock version 2.2");

    mutex = CreateSemaphore(0, 1, 1, 0);
    thread = CreateThread(0, 0, callServerThread,
                          this, 0, &tid);
    if (!thread) THROW("Error starting GUI server thread.");

  } catch (const Exception &e) {
    WSACleanup();
    throw e;
  }
#else
  pthread_mutex_init(&mutex, 0);

  if (pthread_create(&thread, 0, callServerThread, (void *)this))
    THROW("Error starting GUI server thread 2.");
#endif
}

void GUIServer::startUpdate() {
  lock();
}

void GUIServer::endUpdate() {
  // rintf("GUI data served\n");
  request = GS_NO_REQUEST;
  unlock();
}

void GUIServer::lock() {
#ifdef _WIN32
  WaitForSingleObject(mutex, INFINITE);
#else
  pthread_mutex_lock(&mutex);
#endif
}

void GUIServer::unlock() {
#ifdef _WIN32
  ReleaseSemaphore(mutex, 1, 0);
#else
  pthread_mutex_unlock(&mutex);
#endif
}
