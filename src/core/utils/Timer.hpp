#ifndef __UTILS_TIMING_HPP
#define __UTILS_TIMING_HPP

#include <string>
#include <unordered_map>
#include <map>

#include <mpi.h>
#include "utils/statistics/RunningAverage.hpp"

/**
 * @brief Time events and keep staistics.
 *
 * Times are wall time and in seconds.
 */
namespace Utils { namespace Timing {    
    class Timer
    {
    public:
      struct Stats {
	Stats() {}
	Stats(double avg, double sig, int n)
	  : m_avg(avg), m_sig(sig), m_n(n)
	{}

	double avg() const {
	  return m_avg;
	}

	double sig() const {
	  return m_sig;
	}

	int n() const {
	  return m_n;
	}
	
	double m_avg, m_sig;
	int m_n;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
	  ar & m_avg;
	  ar & m_sig;
	  ar & m_n;
	}
      };
      
      void start() {
	m_mark = MPI_Wtime();
      }
      
      double stop() {
	const double time = MPI_Wtime() - m_mark;
	m_running_average.add_sample(time);

	return time;
      }
      
      Stats stats() const {
	return Stats(m_running_average.avg(),
		     m_running_average.sig(),
		     m_running_average.n());
      }

      void reset() {
	m_running_average.clear();
      }
      
      static Timer& get_timer(const std::string &name) {
	return m_timers[name];
      }
      
      static std::map<std::string, Stats> get_stats() {
	std::map<std::string, Stats> ret;
	
	for (auto &it : m_timers) {
	  ret[it.first] = it.second.stats();
	}

	return ret;
      }
      
    private:
      static std::unordered_map<std::string, Timer> m_timers;
      Statistics::RunningAverage<double> m_running_average;
      double m_mark;
    };
  }
}

#endif
