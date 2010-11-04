#ifndef HEMELB_VIS_STREAKLINEDRAWER_H
#define HEMELB_VIS_STREAKLINEDRAWER_H

#include <vector>

#include "constants.h"
#include "mpiInclude.h"
#include "net.h"

#include "lb/GlobalLatticeData.h"
#include "lb/LocalLatticeData.h"

namespace hemelb
{
  namespace vis
  {

    // Class that controls the drawing of streaklines - lines that trace
    // the path of an imaginary particle were it dropped into the fluid.
    class StreaklineDrawer
    {
      public:
        // Constructor and destructor.
        StreaklineDrawer(Net *net, lb::GlobalLatticeData &iGlobLatDat);
        ~StreaklineDrawer();

        // Method to reset streakline drawer
        void Restart();

        // Drawing methods.
        void StreakLines(int time_steps,
                         int time_steps_per_cycle,
                         lb::GlobalLatticeData &iGlobLatDat,
                         lb::LocalLatticeData &iLocalLatDat,
                         Net *net);
        void render(lb::GlobalLatticeData &iGlobLatDat);

      private:

        // Struct for a particle dropped into the fluid.
        struct Particle
        {
            float x, y, z;
            float vx, vy, vz;
            float vel;
            int inlet_id;
        };
        typedef std::vector<Particle> ParticleVector;

        // Struct for information about the velocity field at some point.
        struct VelSiteData
        {
            int proc_id, counter, site_id;
            float vx, vy, vz;
        };

        // Struct to contain the whole velocity field.
        struct VelocityField
        {
            VelSiteData *vel_site_data;
        };

        // Struct for one processor to hold information about its neighbouring processors.
        typedef std::vector<float> FloatVector;
        struct NeighProc
        {
            int id;
            int send_ps, recv_ps;
            int send_vs, recv_vs;

            FloatVector p_to_send, p_to_recv;
            float *v_to_send, *v_to_recv;

            short int *s_to_send, *s_to_recv;
        };

        // Necessary to keep a local store of the number of blocks created, so that we can
        // write a correct constructor. Alternative (TODO) is to use a vector.
        int num_blocks;

        // Counter keeps track of the number of VelSiteDatas created
        int counter;

        // Variables for tracking the actual numbers of particles, and the maximum number
        //(i.e. the number for which memory has been allocated).
        unsigned int nParticles;
        unsigned int nParticleSeeds;
        int particles_to_send_max, particles_to_recv_max;

        // Variables for counting the processors involved etc.
        int neigh_procs;
        int shared_vs;
        int procs;

        // Pointers to the structs.
        VelocityField *velocity_field;
        ParticleVector particleVec;
        ParticleVector particleSeedVec;

        // Arrays for communicating between processors.
        float *v_to_send, *v_to_recv;
        short int *s_to_send, *s_to_recv;
        short int *from_proc_id_to_neigh_proc_index;

        NeighProc neigh_proc[NEIGHBOUR_PROCS_MAX];

        // If using MPI, require these for inter-processor comms.
#ifndef NOMPI
        MPI_Status status[4];
        MPI_Request *req;
#endif

        // Private functions for the creation / deletion of particles.
        void createSeedParticles();
        void createParticle(float x, float y, float z, float vel, int inlet_id);
        void deleteParticle(unsigned int p_index);

        // Private functions for initialising the velocity field.
        void initializeVelFieldBlock(lb::GlobalLatticeData &iGlobLatDat,
                                     int site_i,
                                     int site_j,
                                     int site_k,
                                     int proc_id);
        VelSiteData *velSiteDataPointer(lb::GlobalLatticeData &iGlobLatDat,
                                        int site_i,
                                        int site_j,
                                        int site_k);
        void particleVelocity(Particle *particle_p,
                              float v[2][2][2][3],
                              float interp_v[3]);
        void localVelField(int p_index,
                           float v[2][2][2][3],
                           int *is_interior,
                           lb::GlobalLatticeData &iGlobLatDat,
                           lb::LocalLatticeData &iLocalLatDat,
                           Net *net);

        // Private functions for updating the velocity field and the particles in it.
        void updateVelField(int stage_id,
                            lb::GlobalLatticeData &iGlobLatDat,
                            lb::LocalLatticeData &iLocalLatDat,
                            Net *net);
        void updateParticles();

        // Private functions for inter-proc communication.
        void communicateSiteIds();
        void communicateVelocities(lb::GlobalLatticeData &iGlobLatDat);
        void communicateParticles(lb::GlobalLatticeData &iGlobLatDat, Net *net);

    };

  }
}

#endif // HEMELB_VIS_STREAKLINEDRAWER_H
