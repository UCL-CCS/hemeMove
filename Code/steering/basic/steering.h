#ifndef HEMELB_STEERING_BASIC_STEERING_H
#define HEMELB_STEERING_BASIC_STEERING_H

#include <pthread.h>
#include <semaphore.h>

#include "lb.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace steering
  {

    extern pthread_mutex_t network_buffer_copy_lock;
    extern pthread_mutex_t LOCK;
    extern pthread_cond_t network_send_frame;

    extern sem_t nrl;
    extern sem_t connected_sem;
    extern sem_t steering_var_lock;

    extern bool is_frame_ready;
    extern bool connected;
    extern bool sending_frame;

    extern int send_array_length;


    extern pthread_mutex_t steer_param_lock;

    extern bool updated_mouse_coords;

    void  *hemeLB_steer (void*);
    void UpdateSteerableParameters (int*, vis::Control*, LBM*);
    void* hemeLB_network(void *ptr);
  }
}


#endif // HEME_STEERING_BASIC_STEERING_H
