/// \file timer.h
/// \brief portable implementation of wall clock and cpu timers

#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#include <boost/predef.h>

#if BOOST_OS_WINDOWS
#ifndef NOMINMAX
#define NOMINMAX // otherwise std::max and std::min is defined
#endif
#include <windows.h>
#else
#include <ctime>
#endif

#include <chrono>
#include <stdexcept>

namespace NiHu
{

/// \brief wall clock time
class wc_time
{
	typedef std::chrono::system_clock clock_t;

public:
	/// \brief the time point type returned by tic
	typedef std::chrono::time_point<clock_t> time_point_t;

	/// \brief returns the current time point
	static time_point_t tic()
	{
		return clock_t::now();
	}

	/// \brief returns the time elapsed since t0
	/// \param [in] t0 the reference time point
	/// \return the time elapsed in secs
	static double toc(time_point_t const &t0)
	{
		std::chrono::duration<double, std::nano> dur = tic() - t0;
		return dur.count() / 1e9;
	}
};


/// \brief CPU time
class cpu_time
{
public:
#if BOOST_OS_WINDOWS
	/// \brief the time point type returned by tic
	typedef double time_point_t;

	/// \brief returns the current time point
	static time_point_t tic()
	{
		FILETIME createTime;
		FILETIME exitTime;
		FILETIME kernelTime;
		FILETIME userTime;
		if (GetProcessTimes(GetCurrentProcess(),
			&createTime, &exitTime, &kernelTime, &userTime) == -1)
			throw std::runtime_error("Could not get CPU time point");

		SYSTEMTIME userSystemTime;
		if (FileTimeToSystemTime(&userTime, &userSystemTime) == -1)
			throw std::runtime_error("Could not get CPU time point");

		return (double)userSystemTime.wHour * 3600.0 +
			(double)userSystemTime.wMinute * 60.0 +
			(double)userSystemTime.wSecond +
			(double)userSystemTime.wMilliseconds / 1000.0;
	}

	/// \brief returns the time elapsed since t0
	/// \param [in] t0 the reference time point
	/// \return the time elapsed in secs
	static double toc(time_point_t const &t0)
	{
		return tic() - t0;
	}

#else

	typedef clock_t time_point_t;

	static time_point_t tic()
	{
		return std::clock();
	}

	static double toc(time_point_t const &t0)
	{
		return double(tic() - t0) / CLOCKS_PER_SEC;
	}

#endif
};

} // end of namespace NiHu

#endif
