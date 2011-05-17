#ifndef HEMELB_STEERING_IMAGESENDCOMPONENT_H
#define HEMELB_STEERING_IMAGESENDCOMPONENT_H

#include "constants.h"
#include "lb/lb.h"
#include "lb/SimulationState.h"
#include "lb/LbmParameters.h"
#include "steering/ClientConnection.h"
#include "steering/basic/SimulationParameters.h"

namespace hemelb
{
  namespace steering
  {
    class ImageSendComponent
    {
      public:
        ImageSendComponent(lb::LBM* lbm,
                           lb::SimulationState* iSimState,
                           vis::Control* iControl,
                           const lb::LbmParameters* iLbmParams,
                           ClientConnection* iClientConnection);
        ~ImageSendComponent();

        void DoWork();

        bool ShouldRenderNewNetworkImage();

        bool isFrameReady;
        bool isConnected;
        int send_array_length;

      private:
        ssize_t SendSuccess(int iSocket, char * data, int length);

        ClientConnection* mClientConnection;
        lb::LBM* mLbm;
        lb::SimulationState* mSimState;
        vis::Control* mVisControl;
        const lb::LbmParameters* mLbmParams;

        char* xdrSendBuffer;
        double lastRender;

        // data per pixel is
        // 1 * int (pixel index)
        // 3 * int (pixel RGB)
        static const int bytes_per_pixel_data = 4 * 4;

        // Sent data:
        // 2 * int (pixelsX, pixelsY)
        // 1 * int (bytes of pixel data)
        // pixel data (variable, up to COLOURED_PIXELS_MAX * bytes_per_pixel_data)
        // SimulationParameters::paramsSizeB (metadata - mouse pressure and stress etc)
        static const int maxSendSize = 2 * 4 + 1 * 4 + vis::Screen::COLOURED_PIXELS_MAX
            * bytes_per_pixel_data + SimulationParameters::paramsSizeB;
    };
  }
}

#endif /* HEMELB_STEERING_IMAGESENDCOMPONENT_H */
