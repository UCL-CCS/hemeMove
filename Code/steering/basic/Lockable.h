#ifndef HEMELB_STEERING_BASIC_LOCKABLE
#define HEMELB_STEERING_BASIC_LOCKABLE

#include <semaphore.h>

namespace hemelb
{
  namespace steering
  {

    /**
     * Wrap a variable with locking through semaphores.
     * Construct with Lockable<type>(value).
     * Get a copy of the value with GetValue().
     * Set value with SetValue(newValue).
     * Hold/release the lock with Lock()/Unlock.
     *
     */
    template<typename T>
    class Lockable
    {
      public:
        Lockable(T value) :
          mValue(value)
        {
          sem_init(&mSem, 0, 1);
        }

        ~Lockable()
        {
          sem_destroy(&mSem);
        }

        T GetValue(void)
        {
          sem_wait(&mSem);
          T copy = mValue;
          sem_post(&mSem);
          return copy;
        }

        void SetValue(T val)
        {
          sem_wait(&mSem);
          mValue = val;
          sem_post(&mSem);
        }

        void Lock(void)
        {
          sem_wait(&mSem);
        }

        void Unlock(void)
        {
          sem_post(&mSem);
        }

      private:
        sem_t mSem;
        T mValue;
    };

  } // namespace steering
} // namespace hemelb
#endif // HEMELB_STEERING_BASIC_LOCKABLE
