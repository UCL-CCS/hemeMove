#include <unistd.h>
#include <cerrno>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <csignal>
#include <sys/stat.h>
#include <cstdio>

#include "debug/Debugger.h"
#include "log/Logger.h"
#include "steering/Network.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace steering
  {
    Network::Network(int steeringSessionId) :
      clientConnection(steeringSessionId)
    {

    }

    /**
     * Receive a bytestream of known length from a socket into a buffer.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns true if we have successfully provided that much data.
     */
    bool Network::recv_all(char *buf, const int length)
    {
      // Get a socket
      int socketToClient = clientConnection.GetWorkingSocket();

      if (socketToClient < 0)
      {
        return false;
      }

      ssize_t bytesGot = 0;

      // If we have some buffered data to receive, include that in our count.
      if (recvBuf.length() > 0)
      {
        bytesGot += recvBuf.length();
      }

      // While some data left to be received...
      while (bytesGot < length)
      {
        // Receive some data (up to the remaining length)
        ssize_t n = recv(socketToClient, buf + bytesGot, length - bytesGot, 0);

        // If there was an error, report it and return.
        if (n <= 0)
        {
          // If there was no data and it wasn't simply that the socket would block,
          // raise an error.
          if (errno != EAGAIN)
          {
            log::Logger::Log<log::Warning, log::Singleton>("Steering component: broken network pipe... (%s)",
                                                           strerror(errno));
            Break(socketToClient);
          }
          else
          {
            // If we'd already received some data, store this in the buffer for the start of the
            // next receive.
            if (bytesGot > 0)
            {
              // The newly-received-byte count is the total received byte count minus the length
              // of the buffer.
              long int numNewBytes = bytesGot - recvBuf.length();
              recvBuf.append(buf + recvBuf.length(), numNewBytes);
            }
          }

          // We didn't fully receive.
          return false;
        }
        else
        {
          bytesGot += n;
        }
      }

      // Successfully received what we needed to. Now use the buffer to fill in the gaps, if
      // we were using the buffer at the front of the received data.
      if (recvBuf.length() > 0)
      {
        // The length of buffer to use is the minimum of the buffer size and the length of
        // the string.
        int copyLength = util::NumericalFunctions::min((int) recvBuf.length(), length);

        memcpy(buf, recvBuf.c_str(), copyLength);

        // Empty that much of the buffer.
        recvBuf.erase(0, copyLength);
      }

      return true;
    }

    void Network::PreReceive()
    {
      // Calling send_all will attempt to flush the send buffer.
      send_all(NULL, 0);
    }

    bool Network::IsConnected()
    {
      return clientConnection.GetWorkingSocket() > 0;
    }

    /**
     * Send all bytes from a buffer of known length over a socket.
     *
     * @param sockid
     * @param buf
     * @param length
     * @return Returns the number of bytes sent or -1 on failure.
     */
    bool Network::send_all(const char *buf, const int length)
    {
      // Get a socket.
      int socketToClient = clientConnection.GetWorkingSocket();

      // If there's no such socket, we don't have a connection. Nothing to do.
      if (socketToClient < 0)
      {
        return false;
      }

      // If we have buffered strings to be sent, send those first.
      if (sendBuf.length() > 0)
      {
        long sent = sendInternal(sendBuf.c_str(), sendBuf.length(), socketToClient);

        // Broken socket?
        if (sent < 0)
        {
          return false;
        }

        // Did sending block? We'll try again next time.
        if (sent < (long) sendBuf.length())
        {
          sendBuf.erase(0, sent);

          // If the sending blocks, just ignore this most recent frame. This way, we stop the
          // memory usage of the steering process spiralling, and keep current the images delivered
          // to the client.
          // What we *would* do is sendBuf.append(buf, length);

          return true;
        }
        // If not, we sent the whole buffer.
        else
        {
          sendBuf.clear();
        }
      }

      // If we sent the whole buffer, try to send the new data.
      long sent_bytes = sendInternal(buf, length, socketToClient);

      // Is the socket broken?
      if (sent_bytes < 0)
      {
        return false;
      }
      // Did the socket block? Still return true, because we'll try again next time.
      else if (sent_bytes < length)
      {
        sendBuf.append(buf + sent_bytes, length - sent_bytes);
      }

      // Otherwise we sent the whole thing.
      return true;
    }

    void Network::Break(int socket)
    {
      clientConnection.ReportBroken(socket);

      recvBuf.clear();
      sendBuf.clear();
    }

    /**
     * Returns -1 for a broken socket, otherwise the number of bytes we could send before blocking.
     *
     * @param data
     * @param length
     * @param socket
     * @return
     */
    long Network::sendInternal(const char* data, long length, int socket)
    {
      long sent_bytes = 0;

      while (sent_bytes < length)
      {
        long n = send(socket, data + sent_bytes, length - sent_bytes, 0);

        if (n <= 0)
        {
          // Distinguish between cases where the pipe fails because it'd block
          // (No problem, we'll try again later) or because the pipe is broken.
          if (errno == EWOULDBLOCK)
          {
            return sent_bytes;
          }
          else
          {
            log::Logger::Log<log::Warning, log::Singleton>("Network send had broken pipe... (%s)",
                                                           strerror(errno));
            Break(socket);

            return -1;
          }
        }
        else
        {
          sent_bytes += n;
        }
      }

      return sent_bytes;
    }

  }
}
