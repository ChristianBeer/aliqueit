
#ifndef __misc_h
#define __misc_h

#include <sstream>
#include <string>
#include <vector>
using namespace std;


//changes the size of file <fname> to <size> by either truncating or expanding
void fchsize( string fname, long size );

bool isblank( char c );

//remove surrounding whitespace
string trim( string s );

//make lowercase
string tolower( string s );

template<class _t>
string tostring( _t n ) {
	stringstream ss;
	ss << n;
	return ss.str();
}

string maketwodigit( int n );

//[Dec 12, 12:34:09] running check...
//returns "[Mon DD YYYY, HH:MM:SS]" string of date <timesecs> (returned from time()) or current time if timesecs == 0
string get_timestamp( time_t timesecs = 0 );

//prints a log msg to the log file, with timestamp
void log_msg( string msg, bool timestamp = true );

//prints a msg to both the screen and the log file
void log_and_print( string msg );

//checks if a string consists of digits only
bool isnumber( string s );

//converts 1000000 to 1e6
string scientify( int in );

//appends the contents of file <src_fname> to file <dest_fname>
//it's not the fastest way, but at least the code's portable...
void append_file( string src_fname, string dest_fname );

//sets process priority from 0 to 4 with 0 being idle and 4 being high.
void set_priority( int priority );

//change current working directory
void cd( string new_dir );

//deletes a file
void delete_file( string file );

char tohex( int n );

//encodes special characters in a url string to %xx hex values
string urlencode( string in );


#endif //__misc_h
