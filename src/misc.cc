
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include "cfg.h"
#include "misc.h"

#if defined(WIN32)
#include <windows.h>
#include <direct.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#else	//linux
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#define _chdir chdir
#define _unlink unlink
#define _getcwd getcwd
#endif

using namespace std;


//changes the size of file <fname> to <size> by either truncating or expanding

void fchsize(string fname, long size) {
#if defined(WIN32)
    int fd;
    if (!(_sopen_s(&fd, fname.c_str(), _O_WRONLY, _SH_DENYNO, _S_IREAD | _S_IWRITE))) {
        _chsize(fd, size);
        _close(fd);
    }
#else	//LINUX
    truncate(fname.c_str(), size);
#endif
}

bool isblank(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

//remove surrounding whitespace

string trim(string s) {
    size_t start, end;
    for (start = 0; s.size() > start && isblank(s[start]); ++start) {
    }
    for (end = s.size(); end > start && isblank(s[end - 1]); --end) {
    }
    return s.substr(start, end - start);
}

//make lowercase

string tolower(string s) {
    for (size_t j = 0; j < s.size(); ++j) {
        s[j] = tolower(s[j]);
    }
    return s;
}

string maketwodigit(int n) {
    return ( n < 10 ? "0" : "") +tostring(n);
}

//[Dec 12, 12:34:09] running check...
//returns "[Mon DD YYYY, HH:MM:SS]" string of date <timesecs> (returned from time()) or current time if timesecs == 0

string get_timestamp(time_t timesecs) {
    const char *months[] = {"", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    struct tm * t;
    time_t tt = timesecs;
    if (!tt) tt = time(0);
    //localtime_s( &t, &tt );
    t = localtime(&tt);
    stringstream ss;
    ++t->tm_mon;

    ss << "[" << months[t->tm_mon] << " " << (t->tm_mday < 10 ? "0" : "") << t->tm_mday << " " << (t->tm_year + 1900) << ", ";
    ss << maketwodigit(t->tm_hour) << ":" << maketwodigit(t->tm_min) << ":" << maketwodigit(t->tm_sec) << "]";
    return ss.str();
}

//prints a log msg to the log file, with timestamp

void log_msg(string msg, bool timestamp) {
    ofstream f(cfg.log_file.c_str(), ios::app);
    if (!f.is_open()) {
        cout << "WARNING: couldn't open log file for writing!" << endl;
        return;
    }
    if (timestamp) f << get_timestamp() << " ";
    f << msg;
    f.close();
}

//prints a msg to both the screen and the log file

void log_and_print(string msg) {
    cout << msg;
    while (msg.size() && msg[0] == '\n') {
        log_msg("\n", false);
        msg = msg.substr(1);
    }
    if (msg.size()) {
        log_msg(msg);
    }
}

//checks if a string consists of digits only

bool isnumber(string s) {
    for (size_t j = 0; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') return false;
    }
    return true;
}

//converts 1000000 to 1e6

string scientify(int in) {
    string n = tostring(in);
    size_t o;
    for (o = n.size(); o > 0 && n[o - 1] == '0'; --o) {
    }
    if (o + 3 > n.size()) return n;
    return n.substr(0, o) + "e" + tostring(n.size() - o);
}

//appends the contents of file <src_fname> to file <dest_fname>
//it's not the fastest way, but at least the code's portable...

void append_file(string src_fname, string dest_fname) {
    ifstream fi(src_fname.c_str());
    if (!fi.is_open()) return;
    ofstream fo(dest_fname.c_str(), ios::app);
    if (fo.is_open()) {
        string line;
        while (getline(fi, line)) {
            fo << line << endl;
        }
        fo.close();
    }
    fi.close();
}

//sets process priority from 0 to 4 with 0 being idle and 4 being high.

void set_priority(int priority) {
#if defined(WIN32)
    unsigned int aprio[] = {IDLE_PRIORITY_CLASS, BELOW_NORMAL_PRIORITY_CLASS, NORMAL_PRIORITY_CLASS,
        ABOVE_NORMAL_PRIORITY_CLASS, HIGH_PRIORITY_CLASS};
    //unsigned int atprio[] = { THREAD_PRIORITY_IDLE, THREAD_PRIORITY_LOWEST, THREAD_PRIORITY_BELOW_NORMAL,
    //	THREAD_PRIORITY_NORMAL, THREAD_PRIORITY_ABOVE_NORMAL, THREAD_PRIORITY_HIGHEST };
    SetPriorityClass(GetCurrentProcess(), aprio[priority]);
    //SetThreadPriority( GetCurrentThread(), atprio[0] );
#else	//linux
    setpriority(PRIO_PROCESS, 0, 20 - 10 * priority);
#endif
}

//change current working directory

void cd(string new_dir) {
    _chdir(new_dir.c_str());
}

//deletes a file

void delete_file(string file) {
    _unlink(file.c_str());
}

char tohex(int n) {
    return "0123456789abcdef"[n % 16];
}

//encodes special characters in a url string to %xx hex values

string urlencode(string in) {
    string s;
    for (size_t j = 0; j < in.size(); ++j) {
        if (isalnum(in[j])) {
            s += in[j];
        } else if (in[j] == ' ') {
            s += '+';
        } else {
            s += '%';
            s += tohex(in[j] / 16);
            s += tohex(in[j] % 16);
        }
    }
    return s;
}
