#include "../../include/Logger.h"

Logger::Logger(){
  m_filename = "iscam_runtime.log";
  if(file_exists(m_filename)){
    remove(m_filename.c_str());
  }
  m_out.open(m_filename.c_str(), std::ios::out | std::ios::app);
  // Write header for logfile
  m_start_time = time(0);
  m_stream<<"ISCAM runtime log - "<<ctime(&m_start_time)<<'\n';
}

bool Logger::file_exists(const std::string& filename){
 struct stat buf;
 if (stat(filename.c_str(), &buf) != -1){
   return true;
 }
 return false;
}

Logger::~Logger(){
  // Insert end time and total runtime.
  m_end_time = time(0);
  m_elapsed_time = difftime(m_end_time, m_start_time);
	long hour = long(m_elapsed_time) / 3600;
	long minute = long(m_elapsed_time) % 3600 / 60;
	long second = (long(m_elapsed_time) % 3600) % 60;

  m_stream<<'\n';
	m_stream<<"Start time: "<<ctime(&m_start_time);
	m_stream<<"Finish time: "<<ctime(&m_end_time);
	m_stream<<"Runtime: "<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<'\n';

  m_out<<m_stream.rdbuf();
  m_out.flush();
}
