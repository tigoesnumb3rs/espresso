#ifndef __UTILS_TIMING_HPP
#define __UTILS_TIMING_HPP

#include <string>
#include <unordered_map>

#include <mpi.h>
#include "utils/statistics/RunningAverage.hpp"

namespace Utils { namespace Timing {    
    class Timer
    {
    public:
      void start() {
	m_mark = MPI_Wtime();
      }
      void stop() {
	const double time = MPI_Wtime() - m_mark;
	m_running_average.add_sample(time);
      }
      
      const Statistics::RunningAverage<double> &stats() const {
	return m_running_average;
      }

      static Timer& get_timer(const std::string &name) {
	return m_timers[name];
      }
      
      static std::unordered_map<std::string, Timer> get_all_timers() {
	return m_timers;
      }
    private:
      static std::unordered_map<std::string, Timer> m_timers;
      Statistics::RunningAverage<double> m_running_average;
      double m_mark;
    };
  }
}


#endif
