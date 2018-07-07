
#ifndef __aliqueit_h
#define __aliqueit_h

#include <sstream>
#include <string>
#include <vector>

#pragma warning(push)
#pragma warning(disable : 4244 4800)
#include <gmpxx.h>
#pragma warning(pop)

using namespace std;


//finds factors in a log file
int find_log_factors(string fname, string input_number, string required_on_line, string prefix, vector<mpz_class> & factors, bool save_factorless_log = true);

//finds factors in a gmp-ecm log file and only saves the relevant parts
int find_log_factors_gmp_ecm(string fname, string input_number, string required_on_line, string prefix, vector<mpz_class> & factors);

//checks if any factor in <factors> has >= <min_digits> digits and if so prints a special little hooray msg. :D
void check_for_neat_factors(vector<mpz_class> & factors, unsigned int min_digits);


#endif //__aliqueit_h
