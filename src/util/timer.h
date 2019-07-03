#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED

#if defined( _WIN32 ) | defined( _WIN64 )
#include <windows.h>
#endif

#include <chrono>
#include <exception>

namespace NiHu
{

class wc_time
{
public:
	typedef std::chrono::system_clock clock_t;
	typedef std::chrono::time_point<clock_t> time_point_t;

	static time_point_t tic()
	{
		return clock_t::now();
	}

	static double toc(time_point_t const &t0)
	{
		std::chrono::duration<double, std::nano> dur = tic() - t0;
		return dur.count() / 1e9;
	}
};


class cpu_time
{
public:
#if defined( _WIN32 ) | defined( _WIN64 )
	typedef double time_point_t;

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
#endif

	static double toc(time_point_t const &t0)
	{
		return tic() - t0;
	}
};

} // end of namespace NiHu

#endif
