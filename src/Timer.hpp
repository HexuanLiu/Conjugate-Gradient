//
// This file is part of the course materials for AMATH483/583 at the University of Washington,
// Spring 2017
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//

#if !defined(TIMER_HPP)
#define TIMER_HPP

#include <chrono>

class Timer {
private:
  typedef std::chrono::time_point<std::chrono::system_clock> time_t;

public:
  Timer() : startTime(), stopTime() {}

  time_t start() { return (startTime = std::chrono::system_clock::now()); }
  time_t stop() { return (stopTime = std::chrono::system_clock::now()); }
  double elapsed() { return std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count(); }

private:
  time_t startTime, stopTime;
};

#endif    // TIMER_HPP
