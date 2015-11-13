#ifndef _LOGGER_H
#define _LOGGER_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/stat.h>

#define LOG Logger::instance()

class Logger{
 public:
  static Logger& instance(){static Logger inst;return inst;}
  // For const values, i.e. hard coded strings
  template <class T>
    Logger &operator<<(const T& thing){m_stream<<thing;return *this;}
  // For non-const values, i.e. admb types
  template <class T>
    Logger &operator<<(T& thing){m_stream<<thing;return *this;}
  // For correct instantiation of std::endl template parameters (so you can use 'endl' instead of '\n'
  Logger &operator<<(std::ostream& (*pf) (std::ostream&)){m_stream<<pf;return *this;}
 private:
  Logger();
  ~Logger();
  std::string m_filename;
  std::stringstream m_stream;
  std::ofstream m_out;
  time_t m_start_time;
  time_t m_end_time;
  double m_elapsed_time;
  bool file_exists(const std::string& filename);
};
#endif
