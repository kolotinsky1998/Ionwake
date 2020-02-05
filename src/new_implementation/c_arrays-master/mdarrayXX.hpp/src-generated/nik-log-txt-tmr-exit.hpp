/*########################################################################################################
#                                                                                                        #
#  Here are utility functions for:                                                                       #
#      * Logging to a file                                                                               #
#      * Text files output                                                                               #
#      * Timer                                                                                           #
#      * Exit                                                                                            #
#                                                                                                        #
#  The 3-Clause BSD License:                                                                             #
#                                                                                                        #
#  Copyright 2017 Nikolay Khabarov, International Institute for Applied Systems Analysis (IIASA).        #
#  Redistribution and use in source and binary forms, with or without modification, are permitted        #
#  provided that the following conditions are met:                                                       #
#                                                                                                        #
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions     #
#     and the following disclaimer.                                                                      #
#                                                                                                        #
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions  #
#     and the following disclaimer in the documentation and/or other materials provided with the         #
#     distribution.                                                                                      #
#                                                                                                        #
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse  #
#     or promote products derived from this software without specific prior written permission.          #
#                                                                                                        #
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR        #
#  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND      #
#  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            #
#  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL     #
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     #
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER    #
#  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT     #
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       #
#                                                                                                        #
########################################################################################################*/

/*########################################################################################################
# Acknowledgement:                                                                                       #
# This work is supported by the Synergy Grant 610028 IMBALANCE-P: Effects of phosphorus limitations      #
# on Life, Earth system and Society (Horizon 2020 European Union Funding for Research and Innovation).   #
########################################################################################################*/

#ifndef _nik_include_logging_timer_exit_def_
#define _nik_include_logging_timer_exit_def_ 

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdarg>
#include <cfloat>
#include <sys/fcntl.h>
#include <sys/types.h>
#include <sys/time.h>

///////////////////////////////////////////////////////////////////////////
// utility classes (used by utility functions below, a user should better use the function than a class)
///////////////////////////////////////////////////////////////////////////

void printflog(const char *fmt, ...); // forward ref for logging output to main log file

//// Text file writer class
class NikTextFileWriter {
	public:
		NikTextFileWriter();
		~NikTextFileWriter();
		void open(const char *filepath, const char *filename); // if not opened will print to screen only
		void open(const char *fullfilename); // if not opened will print to screen only
		void openfnfs(const char *fmt, ...); // file name is formatted string
		void printf(const char *fmt, ...);
		void vprintf(const char *fmt, va_list args);
	private:
		FILE *fh;
};

//// Logger class
class NikLogger {
	public:
		NikLogger();
		~NikLogger();
		void open(const char *logpath, const char *logfname); // if not opened will print to screen only
		void open(const char *fnamelog); // if not opened will print to screen only
		void printf(const char *fmt, ...);
		void vprintf(const char *fmt, va_list args);
	private:
		FILE *fh;
};

//// Timer class (using some ideas from http://berenger.eu/blog/2010/09/01/c-time-clock-time-manager-cross-platform-windows-posix/)
class NikTimer {
	public:
		NikTimer(NikLogger *logger); 
		void start();
		void show();
	private:		
		double time_start, time_curr, time_last;
		NikLogger *logger;
		double nik_clock(); // compatible with unix and windows
};

///////////////////////////////////////////////////////////////////////////
// utility functions (to be used directly by a user)
///////////////////////////////////////////////////////////////////////////

void starttimer();

void showtimer();

void openlog(const char *file_out_path, const char *log_file_name);

void openlog(const char *log_file_name_with_path);

void printflog(const char *fmt, ...);

void vprintflog(const char *fmt, va_list args);

void err_exit(const char *fmt, ...);

void dbg_exit(const char *fmt, ...);

#endif
