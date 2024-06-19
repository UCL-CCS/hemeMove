// -*- mode: c++; -*-
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HLBGMYTOOL_OCT_PARALLELAPPLY_H
#define HLBGMYTOOL_OCT_PARALLELAPPLY_H

#include <functional>
#include <future>
#include <thread>
#include <vector>

namespace hemelb::gmytool::oct {

template <class ReturnType, class ArgType>
class ParallelApply {
 public:
  typedef std::function<ReturnType(ArgType)> FunctionType;
  typedef std::future<ReturnType> Future;

  ParallelApply(FunctionType f, int NT = 0, int MQ = 0)
      : func(f),
        nThreads(NT ? NT : std::thread::hardware_concurrency()),
        maxQueue(MQ ? MQ : 4 * NT),
        all_work_queued(false) {
    int np = nThreads;
    workers.reserve(np);
    while (np-- > 0)
      AddWorker();
  }

  Future operator()(ArgType arg) {
    assert(!all_work_queued.load());

    QueueItem qi = std::make_pair(Promise(), arg);
    auto future_ans = qi.first.get_future();
    // enqueue
    while (queue.size() > maxQueue) {
      // sleep
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    LockT lock(queue_access);
    queue.push_back(std::move(qi));
    queue_flag.notify_one();
    return future_ans;
  }

  void Done() { all_work_queued = true; }
  ~ParallelApply() {
    Done();
    for (auto& thr : workers)
      thr.join();
  }

 private:
  typedef std::promise<ReturnType> Promise;
  // typedef std::tuple<ArgTypes...> ArgsTuple;

  typedef std::pair<Promise, ArgType> QueueItem;

  typedef std::unique_lock<std::mutex> LockT;

  FunctionType func;
  int nThreads;
  unsigned maxQueue;
  std::vector<std::thread> workers;

  std::atomic<bool> all_work_queued;

  std::deque<QueueItem> queue;
  std::mutex queue_access;
  std::condition_variable queue_flag;

  void AddWorker() {
    workers.emplace_back([this]() {
      // Thread main
      bool running = true;
      while (running) {
        Promise p;
        ArgType arg;
        {
          // Get work from the queue if it's there
          LockT lock(queue_access);
          if (queue.empty()) {
            if (all_work_queued.load()) {
              running = false;
              continue;
            }
            queue_flag.wait_for(lock, std::chrono::milliseconds(1));
            continue;
          }
          // We have the lock and work in the queue
          p = std::move(queue.front().first);
          arg = queue.front().second;
          queue.pop_front();
        }  // this releases the lock

        // Do the work we just got
        p.set_value(func(arg));
      }

      // if not, are we done (if yes, quit) else sleep
      // if there is work, do it
    });
  }
};

}  // namespace hemelb::gmytool::oct
#endif  // HLBGMYTOOL_OCT_PARALLELAPPLY_H
