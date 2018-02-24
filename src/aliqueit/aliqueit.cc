
//aliqueit, computes aliquot sequences
//Mikael Klasson 2009-2011
//mklasson at googles big mail thing
//http://mklasson.com
//
//Additional code by Greebley, bsquared, and bchaffin.


#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "cfg.h"
#include "misc.h"
#include "aliqueit.h"

#pragma warning(push)
#pragma warning(disable : 4244 4800)
#include <gmpxx.h>
#pragma warning(pop)

using namespace std;


const string version = "v1.12";


vector<unsigned int> trial_primes; //precalced primes for trial factoring

map<mpz_class, unsigned int> seq_values; //seq_value -> index. cached sequence values. Used to detect cycles

bool skip_ecm = false; //cmdline arg "-e" skips ecm on the current sequence iteration
bool use_yafu_factor = false; //cmdline arg "-y" uses yafu for everything after initial testing with trial division and rho

bool factor(mpz_class n, vector<pair<mpz_class, int> > & factors, vector<mpz_class> & new_factors, bool clear_factors = true, int max_cofactor = 0, int max_ecm_level = 0);


//checks if <value> restarts a cycle.
//if it does, the cycle is printed and the program exits. otherwise the value is added for future cycle detection.

void add_and_check_cycle(unsigned int index, mpz_class & value) {
    if (seq_values.count(value)) {
        log_and_print("Cycle detected: index " + tostring(seq_values[value]) + " to " + tostring(index - 1) + ".\n");
        exit(0);
    }
    seq_values[value] = index;
}

//checks if <value> < <seq>, i.e. the sequence merges with an earlier sequence number.

void check_merge(mpz_class & seq, mpz_class & value, int index) {
    if (cfg.detect_merge && value < seq) {
        log_and_print("Sequence merges with earlier sequence " + value.get_str() + " at index " + tostring(index) + ".\n");
        exit(0);
    }
}

void precalc_trial_primes() {
    cout << "Precalcing primes for trial factoring..." << endl;
    mpz_class p(2);
    while (p < cfg.trial_cutoff) {
        trial_primes.push_back(p.get_ui());
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    }
}

//sorts factors and merges any <p,x>,<p,y> into <p,x+y>

void merge_factors(vector<pair<mpz_class, int> > & factors) {
    sort(factors.begin(), factors.end());
    for (size_t j = 1; j < factors.size();) {
        if (factors[j].first == factors[j - 1].first) {
            factors[j - 1].second += factors[j].second;
            factors.erase(factors.begin() + j);
        } else {
            ++j;
        }
    }
}

//finds factors in a log file

int find_log_factors(string fname, string input_number, string required_on_line, string prefix, vector<mpz_class> & factors, bool save_factorless_log) {
    string line, factor;
    int found_factors = 0;

    ifstream ftmp(fname.c_str());
    while (getline(ftmp, line)) {
        size_t pos = line.find(input_number);
        if (pos != line.npos) {
            mpz_class input(input_number, 10);
            while (getline(ftmp, line)) {
                pos = line.find(prefix);
                if (pos != line.npos && line.find(required_on_line) != line.npos) {
                    stringstream ss(line.substr(pos + prefix.size()));
                    ss >> ws >> factor;
                    mpz_class fac(factor, 10);
                    if (mpz_divisible_p(input.get_mpz_t(), fac.get_mpz_t())) {
                        factors.push_back(fac);
                        ++found_factors;
                    } else {
                        log_and_print("WARNING: bad factor found: " + factor + "\n");
                    }
                }
            }
            break;
        }
    }
    ftmp.close();
    if (cfg.save_log && (found_factors || save_factorless_log)) {
        append_file(fname, cfg.log_file);
    }
    return found_factors;
}

//finds factors in a gmp-ecm log file and only saves the relevant parts

int find_log_factors_gmp_ecm(string fname, string input_number, string required_on_line, string prefix, vector<mpz_class> & factors) {
    string line, factor;
    int found_factors = 0;
    string save_to_log;
    int curve = 0;

    ifstream ftmp(fname.c_str());
    while (getline(ftmp, line)) {
        if (line.substr(0, 5) == "Using") {
            save_to_log = ""; //reset what we're going to log
            ++curve;
        }
        save_to_log += line + "\n";
        size_t pos = line.find(prefix);
        if (pos != line.npos && line.find(required_on_line) != line.npos) {
            stringstream ss(line.substr(pos + prefix.size()));
            ss >> ws >> factor;
            factors.push_back(mpz_class(factor, 10));
            ++found_factors;
            save_to_log = "Curve " + tostring(curve) + ":\n" + save_to_log;
            log_msg(save_to_log, false);
        }
    }
    ftmp.close();
    return found_factors;
}

//checks if any factor in <factors> has >= <min_digits> digits and if so prints a special little hooray msg. :D

void check_for_neat_factors(vector<mpz_class> & factors, unsigned int min_digits) {
    for (size_t j = 0; j < factors.size(); ++j) {
        string f = factors[j].get_str();
        if (f.size() >= min_digits && mpz_probab_prime_p(factors[j].get_mpz_t(), 25)) {
            string msg = "*** Neat " + tostring(f.size()) + "-digit factor found: " + f + "\n";
            log_and_print(msg);
        }
    }
}

//returns an estimation of the appropriate factor-digit level to do ecm to for the given input size

float get_ecm_level_gnfs(unsigned int input_digits) {
    //magic formula making ecm_time/nfs_time = 1/4 fit decently for the range c100..c120 on my machine, a core i7.
    //this is equivalent to doing ecm to a factor size level of 32.9 for c100 and 37.6 for c120.
    //32.9 is interpreted to mean up to and including 30 digits plus 2.9/5 of the 35-digit level.
    //return 9.4f + 0.235f * input_digits;
    return cfg.gnfs_k * input_digits + cfg.gnfs_m;
}

//returns an estimation of the appropriate factor-digit level to do ecm to for the given input size

float get_ecm_level_qs(unsigned int input_digits) {
    //magic formula making ecm_time/qs_time = 1/4 fit decently for the range c70..c95 on my machine, a core i7..
    //this is equivalent to doing ecm to a factor size level of 20.1 for c70 and 31.3 for c95.
    //return 0.448f * input_digits - 11.26f;
    return cfg.qs_k * input_digits + cfg.qs_m;
}

//returns an estimation of the appropriate factor-digit level to do ecm to for the given input size
//calls the appropriate *_qs or *_gnfs

float get_ecm_level(unsigned int input_digits) {
    if (input_digits < cfg.gnfs_cutoff) {
        return get_ecm_level_qs(input_digits);
    } else {
        return get_ecm_level_gnfs(input_digits);
    }
}

//run ecm on big input_numbers using ~ gmp-ecm recommended settings and trying to make the ecm time take a constant fraction of the qs/gnfs time.

int run_ecm_big(string input_number, vector<mpz_class> & new_factors, int max_ecm_level) {
    int ecm_settings[][3] = {
        { 20, 11000, 74},
        { 25, 50000, 214},
        { 30, 250000, 430},
        { 35, 1000000, 904},
        { 40, 3000000, 2350},
        { 45, 11000000, 4480},
        { 50, 43000000, 7553},
        { 55, 110000000, 17769},
        { 60, 260000000, 42017},
        { 65, 850000000, 69408}, //quick, someone find a way for this to be useful!
    };
    float max_factor_size = get_ecm_level((unsigned int) input_number.size());

    int prev_level_digits = 0;
    bool done = false;
    for (int level = 0; !done && level < sizeof ( ecm_settings) / (3 * sizeof ( int)); ++level) {
        int b1 = ecm_settings[level][1];
        int curves = ecm_settings[level][2];
        //Added for the -m option to exit ecm early
        if (max_ecm_level > 0 && level >= max_ecm_level) {
            break;
        }
        if (max_factor_size < ecm_settings[level][0]) {
            curves = (int) (curves * (max_factor_size - prev_level_digits) / (ecm_settings[level][0] - prev_level_digits));
            if (!curves) break; //if curves is 0 we'll ecm forever (or until we find a factor, whichever comes first)
            done = true;
        }

        //P-1
        int pb1 = (int) (b1 * 10 * cfg.b1scale_pm1);
        string msg = "c" + tostring(input_number.size()) + ": running P-1 at B1=" + scientify(pb1) + "...";
        log_msg("\n", false);
        log_msg(msg + "\n");
        cout << msg << "              \r" << flush;
        system(("echo " + input_number + " | " + cfg.ecm_cmd + " -pm1 -B2scale "
                + tostring(cfg.b2scale_pm1) + " " + tostring(pb1) + " > " + cfg.ecm_tempfile).c_str());
        int num_facs = find_log_factors_gmp_ecm(cfg.ecm_tempfile, input_number, "Factor found", ": ", new_factors);
        if (num_facs) {
            check_for_neat_factors(new_factors, cfg.neat_factor_limit_pm1);
            return num_facs;
        }

        //P+1
        pb1 = (int) (b1 * 5 * cfg.b1scale_pp1);
        msg = "c" + tostring(input_number.size()) + ": running P+1 x3 at B1=" + scientify(pb1) + "...";
        log_msg("\n", false);
        log_msg(msg + "\n");
        cout << msg << "              \r" << flush;
        system(("echo " + input_number + " | " + cfg.ecm_cmd + " -one -pp1 -c 3 -B2scale "
                + tostring(cfg.b2scale_pp1) + " " + tostring(pb1) + " > " + cfg.ecm_tempfile).c_str());
        num_facs = find_log_factors_gmp_ecm(cfg.ecm_tempfile, input_number, "Factor found", ": ", new_factors);
        if (num_facs) {
            check_for_neat_factors(new_factors, cfg.neat_factor_limit_pp1);
            return num_facs;
        }

        //ECM
        b1 = (int) (b1 * cfg.b1scale_ecm);
        msg = "c" + tostring(input_number.size()) + ": running " + tostring(curves) + " ecm curves at B1=" + scientify(b1) + "...";
        log_msg("\n", false);
        log_msg(msg + "\n");
        cout << msg << "              \r" << flush;
        if (cfg.use_ecmpy) { //use ecm.py for multithreading
            system(("echo " + input_number + " | " + cfg.ecmpy_cmd + " -one -c " + tostring(curves)
                    + " -B2scale " + tostring(cfg.b2scale_ecm) + " -out " + cfg.ecm_tempfile + " " + tostring(b1)).c_str());
        } else { //use regular ecm
            system(("echo " + input_number + " | " + cfg.ecm_cmd + " -one -c " + tostring(curves)
                    + " -B2scale " + tostring(cfg.b2scale_ecm) + " " + tostring(b1) + " > " + cfg.ecm_tempfile).c_str());
        }
        num_facs = find_log_factors_gmp_ecm(cfg.ecm_tempfile, input_number, "Factor found", ": ", new_factors);
        if (num_facs) {
            check_for_neat_factors(new_factors, cfg.neat_factor_limit_ecm);
            return num_facs;
        }

        prev_level_digits = ecm_settings[level][0];
    }
    return 0;
}

//run auto-increasing ecm using #curves decided by <ecm_depth> or infinitely many if <loop_forever> is true

int run_ecm_autoinc(string input_number, vector<mpz_class> & new_factors, bool loop_forever = false) {
    delete_file(cfg.ecm_tempfile);

    //find #curves to run
    size_t input_digits = input_number.size();
    size_t digits = 0, curves = 0;
    if (!loop_forever) {
        string line;
        stringstream ss(cfg.ecm_depth);
        while (ss >> digits) {
            if (digits > input_digits) break;
            getline(ss, line, ',');
            curves = atoi(line.c_str());
        }

        if (!curves) return 0;
    }

    string msg = "c" + tostring(input_digits) + ": running " + tostring(curves) + " auto-increasing ecm curves...";
    log_msg("\n", false);
    log_msg(msg + "\n");
    cout << msg << "                 \r" << flush;

    system(("echo " + input_number + " | " + cfg.ecm_cmd + " -one -c " + tostring(curves) + " -I 1 1 > " + cfg.ecm_tempfile).c_str());

    int num_facs = find_log_factors_gmp_ecm(cfg.ecm_tempfile, input_number, "Factor found", ": ", new_factors);
    check_for_neat_factors(new_factors, cfg.neat_factor_limit_ecm);
    return num_facs;
}

int run_ecm(string input_number, vector<mpz_class> & new_factors, int max_ecm_level) {
    if (input_number.size() < cfg.big_ecm_cutoff) {
        return run_ecm_autoinc(input_number, new_factors);
    } else {
        return run_ecm_big(input_number, new_factors, max_ecm_level);
    }
}

//do at most <max_iterations> steps of Pollard's rho algorithm using f(x)=x^2+<c> (mod <n>)
//returns #factors found (if <new_factors> was empty)

int rho(mpz_class & n, vector<mpz_class> & new_factors, unsigned int c, int max_iterations) {
    mpz_t x, y, d, t;
    mpz_init_set_ui(x, 2);
    mpz_init_set_ui(y, 2);
    mpz_init_set_ui(d, 1);
    mpz_init_set_ui(t, 1);
    mpz_ptr pn = n.get_mpz_t();
    unsigned int pow = 1, lam = 1;
    for (int j = 1; j <= max_iterations; ++j) {
        //Floyd's cycle finder
        //mpz_powm_ui( x, x, 2, pn );
        //mpz_add_ui( x, x, c );
        //mpz_powm_ui( y, y, 2, pn );
        //mpz_add_ui( y, y, c );
        //mpz_powm_ui( y, y, 2, pn );
        //mpz_add_ui( y, y, c );
        //mpz_sub( d, x, y );
        //mpz_abs( d, d );

        //Brent's cycle finder
        if (pow == lam) {
            mpz_set(x, y);
            pow *= 2;
            lam = 0;
        }
        ++lam;
        mpz_powm_ui(y, y, 2, pn);
        mpz_add_ui(y, y, c);
        mpz_sub(d, x, y);
        mpz_abs(d, d);

        //gcd phase
        mpz_mul(t, t, d);
        mpz_mod(t, t, pn);
        if (j % 100 == 0 || j == max_iterations) {
            mpz_gcd(t, t, pn);
            if (mpz_cmp_ui(t, 1)) {
                if (mpz_cmp(t, pn)) {
                    new_factors.push_back(mpz_class(t));
                }
                break;
            }
        }
    }
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(d);
    mpz_clear(t);
    return (int) new_factors.size();
}

//returns #factors found (only finds max 1 factor 1 now)

int run_rho(mpz_class & input_number, vector<mpz_class> & new_factors) {
    log_msg("\n", false);
    log_msg("c" + tostring(input_number.get_str().size()) + ": running rho...\n");

    //use yafu's rho
    //delete_file( "factor.log" );
    //system( ( yafu_cmd + " \"rho(" + input_number + ")\" > " + null_device ).c_str() );
    //return find_log_factors( "factor.log", input_number, ": found ", " = ", new_factors );

    //use internal rho for significant speed increase on small numbers (no log parsing etc)
    if (rho(input_number, new_factors, 1, 15000)) {
        return 1;
    }
    //if( rho( input_number, new_factors, 3, 6000 ) ) {
    //	return 1;
    //}
    //if( rho( input_number, new_factors, 2, 6000 ) ) {
    //	return 1;
    //}
    return 0;
}

//returns #factors found

int run_yafu(string input_number, vector<mpz_class> & new_factors) {
    delete_file("factor.log"); //remove old log...

    string hide_output = "";
    if (input_number.size() < cfg.hide_output_limit) {
        hide_output = " > " + cfg.null_device; //avoid spamming the screen when factorisation just takes a couple of seconds
        cout << "c" << input_number.size() << ": running qs (yafu)...                                         \r" << flush;
    }

    log_msg("\n", false);
    log_msg("c" + tostring(input_number.size()) + ": running qs (yafu)...\n");

    system((cfg.yafu_cmd + " \"siqs(" + input_number + ")\"" + hide_output).c_str());

    int num_found = find_log_factors("factor.log", input_number, " prp", " = ", new_factors);
    if (!num_found) {
        cout << "WARNING: couldn't find factor in factor.log" << endl;
    }
    return num_found;
}

//returns #factors found

int run_factor(string input_number, vector<mpz_class> & new_factors, bool no_ecm, bool pretest_only) {
    delete_file("factor.log"); //remove old log...

    string hide_output = "";
    if (input_number.size() < cfg.hide_output_limit) {
        hide_output = " > " + cfg.null_device; //avoid spamming the screen when factorisation just takes a couple of seconds
        cout << "c" << input_number.size() << ": running factor (yafu)...                                         \r" << flush;
    }

    log_msg("\n", false);
    log_msg("c" + tostring(input_number.size()) + ": running factor (yafu)...\n");

    string pretest_arg = pretest_only ? " -pretest " : "";
    string noecm_arg = no_ecm ? " -noecm " : "";

    system((cfg.yafu_cmd + " \"factor(" + input_number + ")\"" + pretest_arg + noecm_arg + hide_output).c_str());

    int num_found = find_log_factors("factor.log", input_number, " prp", " = ", new_factors);
    if (!num_found && !pretest_only) {
        cout << "WARNING: couldn't find factor in factor.log" << endl;
    }
    return num_found;
}


//returns #factors found

int run_msieve(string input_number, vector<mpz_class> & new_factors) {
    delete_file("msieve.log"); //remove old log...

    string hide_output = "";
    if (input_number.size() < cfg.hide_output_limit) {
        hide_output = " 2> " + cfg.null_device; //avoid spamming the screen when factorisation just takes a couple of seconds
        cout << "c" << input_number.size() << ": running qs (msieve)...                                         \r" << flush;
    }

    log_msg("\n", false);
    log_msg("c" + tostring(input_number.size()) + ": running qs (msieve)...\n");

    system((cfg.msieve_cmd + " " + input_number + hide_output).c_str());
    int num_found = find_log_factors("msieve.log", input_number, "factor:", "factor:", new_factors);
    if (!num_found) {
        cout << "WARNING: couldn't find factor in msieve.log" << endl;
    }
    return num_found;
}

//run qs using yafu/msieve (or both if the first fails for some reason)
//returns #factors found

int run_qs(string input_number, vector<mpz_class> & new_factors) {
    int num_found;
    if (cfg.prefer_yafu) {
        if (!(num_found = run_yafu(input_number, new_factors))) {
            num_found = run_msieve(input_number, new_factors);
        }
    } else {
        if (!(num_found = run_msieve(input_number, new_factors))) {
            num_found = run_yafu(input_number, new_factors);
        }
    }
    return num_found;
}

//converts from msieve's .fb to ggnfs' .poly format
//returns false if some error occurred

bool convert_poly(string infile, string poly_name) {
    ifstream fi(infile.c_str());
    if (!fi.is_open()) return false;
    string outfile = poly_name + ".poly";
    ofstream fo(outfile.c_str());
    if (!fo.is_open()) {
        fi.close();
        return false;
    }

    fo << "name: " << poly_name << endl;
    fo << "type: gnfs" << endl;

    map<string, pair<string, bool> > args; //map<msieve_arg_name, pair<ggnfs_arg_name, arg_is_ok> >
    args["N"] = make_pair("n", false);
    args["SKEW"] = make_pair("skew", false);
    args["R0"] = make_pair("Y0", false);
    args["R1"] = make_pair("Y1", false);
    args["A0"] = make_pair("c0", false);
    args["A1"] = make_pair("c1", true);
    args["A2"] = make_pair("c2", true);
    args["A3"] = make_pair("c3", true);
    args["A4"] = make_pair("c4", true);
    args["A5"] = make_pair("c5", true); //A5 is optional (quartics are ok)
    args["A6"] = make_pair("c6", true); //A6 is optional (sextics are ok)

    string line, name;
    while (fi >> name) {
        getline(fi, line);
        if (args.count(name)) {
            fo << args[name].first << ": " << line << endl;
            args[name].second = true; //this arg is ok
        } else {
            fo << "# " << name << " " << line << endl; //unknown line. copy it as a comment
        }
    }

    fi.close();
    fo.close();

    for (map < string, pair < string, bool> >::iterator i = args.begin(); i != args.end(); ++i) {
        if (!i->second.second) {
            cout << "WARNING: didn't find required parameter: " << i->first << endl;
            delete_file(outfile); //delete this poly file and let ggnfs find its own poly
            return false;
        }
    }

    return true;
}

//runs msieve to find a gnfs poly

void run_msieve_polyfind(string input_number, string poly_name) {
    system((cfg.msieve_cmd + " -v -np " + input_number).c_str());
    convert_poly("msieve.fb", poly_name);
}

//run gnfs using ggnfs, optionally first running msieve to find a poly

int run_gnfs(string input_number, vector<mpz_class> & new_factors) {
    log_msg("\n", false);
    log_msg("c" + tostring(input_number.size()) + ": running gnfs (ggnfs)...\n");

    string dir = "ggnfs_" + input_number;
    system((cfg.makedir_cmd + " " + dir).c_str());
    cd(dir);
    ofstream fo("test.n");
    fo << "n: " << input_number << endl;
    fo.close();

    if (cfg.use_msieve_polyfind) {
        run_msieve_polyfind(input_number, "test");
    }

    system((cfg.ggnfs_cmd + " " + dir + "/test").c_str());
    int num_found = find_log_factors("ggnfs.log", input_number, "factor:", "factor:", new_factors); //factMsieve.pl syntax
    if (!num_found) {
        num_found = find_log_factors("test.log", input_number, "factor:", "factor:", new_factors); //factmsieve.py syntax
    }
    if (!num_found) {
        num_found = find_log_factors("ggnfs.log", input_number, "-> p: ", "-> p: ", new_factors); //factLat.pl syntax
    }
    if (num_found) {
        system(cfg.ggnfs_clean_cmd.c_str()); //remove work files
    }
    cd("..");
    return num_found;
}

//prints a log msg and adds <factor> to <factors>

void found_factor(mpz_class & factor, vector<pair<mpz_class, int> > & factors) {
    log_msg("*** prp" + tostring(factor.get_str().size()) + " = " + factor.get_str() + "\n");
    factors.push_back(make_pair(factor, 1));
}

void log_cofactor(mpz_class & n) {
    if (n != 1) {
        log_msg("\n", false);
        log_msg("Cofactor " + n.get_str() + " (" + tostring(n.get_str().size()) + " digits)\n");
    }
}

//tests <new_factors> against <n> and if valid adds them to <factors> and removes them from <n>

void add_factors(mpz_class & n, vector<pair<mpz_class, int> > & factors, vector<mpz_class> & new_factors) {
    vector<mpz_class> empty_dummy;
    for (size_t j = 0; j < new_factors.size(); ++j) {
        if (!mpz_divisible_p(n.get_mpz_t(), new_factors[j].get_mpz_t())) {
            //check if gcd>1
            mpz_gcd(new_factors[j].get_mpz_t(), n.get_mpz_t(), new_factors[j].get_mpz_t());
            if (new_factors[j] == 1) {
                log_and_print("WARNING: factor doesn't divide n\n");
                continue;
            }
        }
        mpz_divexact(n.get_mpz_t(), n.get_mpz_t(), new_factors[j].get_mpz_t());
        if (mpz_probab_prime_p(new_factors[j].get_mpz_t(), 25)) {
            found_factor(new_factors[j], factors);
        } else {
            log_msg("*** c" + tostring(new_factors[j].get_str().size()) + " = " + new_factors[j].get_str() + "\n");
            factor(new_factors[j], factors, empty_dummy, false); //factor the found factor further and add its prime components to <factors>
        }
    }
    new_factors.clear();
    log_cofactor(n);
}

//factors <n> and returns its prime components (with exponents) in <factors>.
//vector<pair<p_i,x_i> >, n = product(p_i^x_i)

bool factor(mpz_class n, vector<pair<mpz_class, int> > & factors, vector<mpz_class> & new_factors, bool clear_factors, int max_cofactor, int max_ecm_level) {
    if (clear_factors) {
        factors.clear();
    }
    if (n == 1 || mpz_probab_prime_p(n.get_mpz_t(), 25)) {
        found_factor(n, factors);
        return true;
    }
    if (new_factors.size()) {
        add_factors(n, factors, new_factors);
    }

    bool trial_success = false;
    for (size_t j = 0; j < trial_primes.size(); ++j) {
        while (mpz_divisible_ui_p(n.get_mpz_t(), trial_primes[j])) {
            mpz_class tp(trial_primes[j]);
            found_factor(tp, factors);
            mpz_divexact_ui(n.get_mpz_t(), n.get_mpz_t(), trial_primes[j]);
            trial_success = true;
        }
    }
    if (trial_success) {
        log_cofactor(n);
    }

    if (use_yafu_factor) {
        while (n != 1 && !mpz_probab_prime_p(n.get_mpz_t(), 25)) {
            new_factors.clear();
            string input = n.get_str();
            if (run_rho(n, new_factors)) {
                //found rho factor
            } else if (max_cofactor > 0 && (int) input.size() >= max_cofactor) {
                if (run_factor(input, new_factors, skip_ecm, true) == 0) {
                    log_msg("*** Max cofactor found... " + input + "\n");
                    return false;
                }
            } else if (run_factor(input, new_factors, skip_ecm, false)) {
                // yafu factored the number
            } else if (cfg.stop_on_failure) {
                cout << "WARNING: yafu failed to find a factor. This really shouldn't happen." << endl;
                log_msg("*** yafu failed to find a factor. Ending.\n");
                return false;
            } else {
                run_ecm_autoinc(input, new_factors, true);
            }
            add_factors(n, factors, new_factors);
        }
    } else {

        while (n != 1 && !mpz_probab_prime_p(n.get_mpz_t(), 25)) {
            new_factors.clear();
            string input = n.get_str();
            if (run_rho(n, new_factors)) {
                //found rho factor
            } else if (!skip_ecm && run_ecm(input, new_factors, max_ecm_level)) {
                //found ecm factor
            } else if (max_cofactor > 0 && (int) input.size() >= max_cofactor) {
                log_msg("*** Max cofactor found... " + input + "\n");
                return false;
            } else if (input.size() < cfg.gnfs_cutoff) {
                if (!run_qs(input, new_factors)) {
                    if (!cfg.stop_on_failure) {
                        cout << "WARNING: qs failed to find a factor. This really shouldn't happen." << endl
                                << "I'll cheerfully loop and try again though..." << endl;
                    } else {
                      cout << "WARNING: qs failed to find a factor. This really shouldn't happen." << endl;
                      log_msg("*** qs failed to find a factor. Ending.\n");
                      return false;
                    }
                }
            } else { //bring out the big gnfs guns then
                if (!run_gnfs(input, new_factors)) {
                    if (!cfg.stop_on_failure) {
                        cout << "WARNING: gnfs failed to find a factor. This really shouldn't happen." << endl
                                << "I'll just run ecm till the end of time or a factor turns up..." << endl
                                << "Let's hope you don't run out of disk space before either of those." << endl;
                        run_ecm_autoinc(input, new_factors, true);
                    } else {
                        cout << "WARNING: gnfs failed to find a factor. This really shouldn't happen." << endl;
                        log_msg("*** gnfs failed to find a factor. Ending.\n");
                        return false;
                    }
                }
            }
            add_factors(n, factors, new_factors);
        }
    }

    if (n != 1) {
        found_factor(n, factors); //remaining cofactor
    }

    merge_factors(factors); //sorts factors and merges any <p,x>,<p,y> into <p,x+y>
    return true;
}

//get the elf filename corresponding to a sequence number

string get_elf_name(mpz_class & sequence) {
    return cfg.result_file_prefix + sequence.get_str() + ".elf";
}

//calculates <n>=product(<factors>) and <s>=sigma(<n>)
//also verifies that all factors are prime (or rather strongly probably prime)

void sigma(vector<pair<mpz_class, int> > & factors, mpz_class & s, mpz_class & n) {
    mpz_class t, tmp;
    n = 1;
    s = 1;
    for (size_t j = 0; j < factors.size(); ++j) {
        if (!mpz_probab_prime_p(factors[j].first.get_mpz_t(), 25)) {
            log_and_print("\nERROR: factor not prime: " + factors[j].first.get_str() + "\n");
            exit(1);
        }
        t = 1;
        for (int k = 1; k <= factors[j].second; ++k) {
            mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), k);
            t += tmp;
        }
        s *= t;

        mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second);
        n *= tmp;
    }
}

//verifies the integrity of an elf file and updates <index> and <n> for continuing at the end of the file

void parse_elf(mpz_class & seq, int & index, mpz_class & n, mpz_class & last_n, vector<mpz_class> & external_factors, bool only_verify) {
    log_msg("\n\n", false);
    log_and_print("Verifying elf " + seq.get_str() + "\n");
    ifstream f(get_elf_name(seq).c_str());
    if (!f.is_open()) {
        log_and_print("WARNING: couldn't open elf file.\n");
        exit(1);
    }

    streampos line_start_offset = 0;
    index = 0;
    mpz_class vn, smn = seq;
    string line, s, str_n;
    int prev_index = 0;
    bool first_line = true;
    vector<pair<mpz_class, int> > factors;
    while (getline(f, line)) {
        stringstream ss(line);
        char c;
        ss >> index >> c >> str_n >> c;
        n.set_str(str_n, 10);

        cout << "\rVerifying index " << index << "... " << flush;

        string error_header = "\nERROR: @index " + tostring(index) + ": ";
        string warning_header = "\nWARNING: @index " + tostring(index) + ": ";

        if (first_line) {
            if (index) {
                log_and_print(warning_header + "elf doesn't begin with index 0. Can't verify correctness of first line.\n");
            } else if (n != seq) {
                log_and_print(error_header + "value doesn't match sequence number.\n");
                exit(1);
            }
        } else {
            if (prev_index + 1 != index) {
                log_and_print(error_header + "index != prev_index + 1\n");
                exit(1);
            }
            if (smn != n) {
                log_and_print(error_header + "value != sigma - n\n");
                exit(1);
            }
        }
        first_line = false;

        factors.clear();
        while (ss >> ws >> s) {
            if (s == "*") continue;
            size_t o;
            if ((o = s.find('^')) != s.npos) {
                string p = s.substr(0, o);
                string e = s.substr(o + 1);
                if (!isnumber(p) || !isnumber(e)) {
                    log_and_print(error_header + "not a number: " + s + "\n");
                    exit(1);
                }
                factors.push_back(make_pair(mpz_class(p), atoi(e.c_str())));
            } else {
                if (!isnumber(s)) {
                    log_and_print(error_header + "not a number: " + s + "\n");
                    exit(1);
                }
                factors.push_back(make_pair(mpz_class(s), 1));
            }
        }
        sigma(factors, smn, vn);
        smn -= vn;

        if (vn != n) {
            if (!getline(f, line)) {
                f.close();
                if (only_verify) {
                    log_and_print(warning_header + "partial last line.\n");
                    if (cfg.verify_terminations) {
                        exit(1);
                    }
                    return; //don't truncate file if we're only verifying
                }
                log_and_print(warning_header + "partial last line. Resuming...\n");
                //add factors to external_factors
                for (size_t j = 0; j < factors.size(); ++j) {
                    for (int k = 0; k < factors[j].second; ++k) {
                        external_factors.push_back(factors[j].first);
                    }
                }
                fchsize(get_elf_name(seq), (long) line_start_offset);
                return;
            } else {
                log_and_print(error_header + "product(factors) != value\n");
                exit(1);
            }
        }

        add_and_check_cycle(index, n);
        check_merge(seq, n, index);

        prev_index = index;
        ++index; //start computation at next line if this is the last
        line_start_offset = f.tellg();
    }
    f.close();
    last_n = n;
    n = smn;
    if (only_verify && cfg.verify_terminations && n != 1) {
        log_and_print("\nERROR: @index " + tostring(index-1) + ": " + "sigma(n) != 1\n");
        exit(1);
    }
    log_and_print("Elf " + seq.get_str() + " verified OK!\n");
}

//output a result line to screen and elf file

void save_result(mpz_class & sequence, string s) {
    ofstream f(get_elf_name(sequence).c_str(), ios::app);
    if (!f.is_open()) {
        cout << "ERROR: couldn't open elf file for writing!" << endl;
        return;
    }
    f << s;
    f.close();
}

//scans through the last 100KB of the log file looking for matching factors for our current number

void find_previous_factors(mpz_class & n, vector<mpz_class> & factors) {
    ifstream f(cfg.log_file.c_str());
    if (!f.is_open()) {
        return;
    }
    f.seekg(-102400, ios_base::end); //parse the last 100KB of the log
    if (f.tellg() == (streampos) - 1) {
        f.seekg(0); //start at the beginning then
    }
    string line;
    bool scanning = false;
    mpz_class factor;
    string sn = n.get_str();
    while (getline(f, line)) {
        if (line.find("*** Starting ") != line.npos) {
            scanning = line.find(" " + sn + " ") != line.npos;
        } else if (scanning && (line.find("*** p") != line.npos || line.find("*** c") != line.npos)) {
            //TODO: if this is a composite factor then only use it if we fail to find its subfactors
            size_t o = line.find(" = ");
            if (o != line.npos) {
                factor.set_str(line.substr(o + 3), 10);
                if (mpz_divisible_p(n.get_mpz_t(), factor.get_mpz_t())) {
                    factors.push_back(factor);
                    cout << "using previously found factor " << factor.get_str() << endl;
                }
            }
        }
    }
    f.close();
}

void print_help() {
    cout << "aliqueit " << version << " by Mikael Klasson 2009-2010, http://mklasson.com" << endl
            << "usage: aliqueit [-i <index> <start_value>] [-f <factor>] [-s <start>] [-y] [-e] [-p] [-t] [-u] [-q]" << endl
            << "                [-d <max_digits>] [-c <max_cofactor>] [-r <add_digits>] [-m <max_ecm_level>] [-b] <seq>" << endl
            << "  Computes the aliquot sequence starting at <seq>." << endl
            << "  <factor> will be tested against the next sequence value." << endl
            << "    Multiple factors can be given with multiple -f." << endl
            << "  If no <index> is given the program continues where the elf file ends." << endl
            << "    In this case the elf file is also checked for correctness at startup." << endl
            << "  -y uses YAFU for all factoring after rho" << endl
            << "  -e skips ecm on the current iteration and goes straight to qs/gnfs." << endl
            << "  -p runs at idle priority." << endl
            << "  -q quits after factoring the first number." << endl
            << "  -d quits after sequence reaches <max_digits> digits." << endl
            << "  -c quits after sequence passes ecm with a cofactor of <max_cofactor> or more digits." << endl
            << "  -b with -c and -d quits after a number meets both limits." << endl
            << "  -r with -c or -d (or both) increases max digit/cofactor by <add_digits> if no driver." << endl
            << "  -m ecm will only run to level <max_ecm_level> where '1'=20-digit factor, '2'=25, etc." << endl
            << "  -s submits the iterations <start>+ from <seq>'s elf file to Syd's DB." << endl
            << "    You'll need wget to use -s." << endl
            << "  -t quits after verifying elf file integrity." << endl
            << "  -u verifies elf file integrity but also checks for termination (overrides verify_terminations from config file)." << endl;
}

//submits the elf file belonging to sequence <seq> to Syd's DB.
//only sequence iterations >= <from_iteration> are sent.

bool submit_elf(mpz_class & seq, int from_iteration) {
    ifstream f(get_elf_name(seq).c_str());
    if (!f.is_open()) {
        cout << "ERROR: couldn't open elf file: " << get_elf_name(seq) << endl;
        return false;
    }
    string line;
    string post_data = "report=true&msub=";
    //string post_data = "report=true&background=on&msub=";
    int sent_lines = 0;
    while (getline(f, line)) {
        if (atoi(line.c_str()) >= from_iteration) {
            post_data += urlencode(line + "\n");
            ++sent_lines;
        }
    }
    f.close();

    if (!sent_lines) {
        cout << "WARNING: no lines to send." << endl;
        return false;
    }

    string post_file = "aliqueit_tmp_post_" + seq.get_str();
    string tmp_file = "aliqueit_tmp_wget_" + seq.get_str();
    ofstream fo(post_file.c_str());
    fo << post_data;
    fo.close();
    cout << "Sending " << sent_lines << " lines..." << endl;
    system(("wget --cache=off --output-document=" + tmp_file + " --post-file=" + post_file + " http://factordb.com/report.php").c_str());
    delete_file(post_file);

    f.clear();
    f.open(tmp_file.c_str());
    while (getline(f, line)) {
        size_t o = line.find("Found ");
        if (o != line.npos) {
            cout << "DB reports " << atoi(line.substr(o + 6).c_str()) << " new factors." << endl;
        }
    }
    f.close();
    delete_file(tmp_file);
    return true;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        print_help();
        return 1;
    }

    cfg.read_config_file();

    vector<mpz_class> external_factors; //given on cmdline or found in the log
    mpz_class seq = 0, n = 0, last_n = 0;
    int index = 0;
    bool submit = false;
    int submit_start = 0;
    bool quit_after_first = false; //quit after factoring the first number?
    bool quit_after_verify = false; //quit after verifying elf file integrity?
    int max_cofactor = 0;
    int max_digits = 0;
    int max_ecm_level = 0;
    bool both_digits_and_cofactor = false;
    bool exit_only_with_driver = false;
    int no_driver_extra = 0;

    for (int j = 1; j < argc;) {
        if (argv[j][0] == '-') {
            if (argv[j][1] == 'e') {
                skip_ecm = true;
                ++j;
            } else if (argv[j][1] == 'h') {
                print_help();
                return 1;
            } else if (argv[j][1] == 'y') {
                use_yafu_factor = true;
                ++j;
            } else if (argv[j][1] == 'p') {
                set_priority(0);
                ++j;
            } else if (argv[j][1] == 'q') {
                quit_after_first = true;
                ++j;
            } else if (argv[j][1] == 't') {
                quit_after_verify = true;
                ++j;
            } else if (argv[j][1] == 'u') {
                quit_after_verify = true;
                cfg.verify_terminations = true;
                ++j;
            } else if (argv[j][1] == 's') { //submit elf to Syd's DB
                submit = true;
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no number given after " << argv[j] << endl;
                    return 1;
                }
                submit_start = atoi(argv[j + 1]);
                j += 2;
            } else if (argv[j][1] == 'f') {
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no number given after " << argv[j] << endl;
                    return 1;
                }
                external_factors.push_back(mpz_class(argv[j + 1]));
                j += 2;
            } else if (argv[j][1] == 'i') {
                if (j + 2 >= argc || !isnumber(argv[j + 1]) || !isnumber(argv[j + 2])) {
                    cout << "ERROR: bad/no numbers given after " << argv[j] << endl;
                    return 1;
                }
                index = atoi(argv[j + 1]);
                n.set_str(argv[j + 2], 10);
                j += 3;
            } else if (argv[j][1] == 'c') {
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no numbers given after " << argv[j] << endl;
                    return 1;
                }
                max_cofactor = atoi(argv[j + 1]);
                j += 2;
            } else if (argv[j][1] == 'd') {
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no numbers given after " << argv[j] << endl;
                    return 1;
                }
                max_digits = atoi(argv[j + 1]);
                j += 2;
            } else if (argv[j][1] == 'm') {
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no numbers given after " << argv[j] << endl;
                    return 1;
                }
                max_ecm_level = atoi(argv[j + 1]);
                j += 2;
            } else if (argv[j][1] == 'b') {
                both_digits_and_cofactor = true;
                j++;
            } else if (argv[j][1] == 'r') {
                if (j + 1 >= argc || !isnumber(argv[j + 1])) {
                    cout << "ERROR: bad/no numbers given after " << argv[j] << endl;
                    return 1;
                }
                exit_only_with_driver = true;
                no_driver_extra = atoi(argv[j + 1]);
                j += 2;
            } else {
                cout << "ERROR: unknown switch: " << argv[j] << endl;
                return 1;
            }
        } else {
            seq.set_str(argv[j], 10);
            ++j;
        }
    }
    if (seq == 0) {
        cout << "ERROR: no <seq> specified." << endl;
        return 1;
    }
    if (submit) {
        submit_elf(seq, submit_start);
        return 0;
    }
    if (index && n == 0) {
        cout << "ERROR: index given, but no start value." << endl;
        return 1;
    }
    if (quit_after_verify) {
        if (index || n != 0) {
            cout << "ERROR: verify mode doesn't accept index parameters." << endl;
            return 1;
        }
        if (cfg.verify_terminations && cfg.detect_merge) {
          cout << "ERROR: verify_terminations mode needs detect_merge to be off." << endl;
          return 1;
        }
        parse_elf(seq, index, n, last_n, external_factors, true); //verify elf file
        return 0; //return OK. Any error will have triggered an exit( 1 ) already
    }
    if (!index && n == 0) {
        n = seq;
        parse_elf(seq, index, n, last_n, external_factors, false); //verify elf and resume from end...
    }
    if (!index && n != seq) {
        cout << "ERROR: index 0 must have start value == <seq>." << endl;
        return 1;
    }

    precalc_trial_primes();

    cout << "seq = " << seq.get_str() << endl;
    cout << "index = " << index << endl;
    cout << "value = " << n.get_str() << " (" << n.get_str().size() << " digits)" << endl;

    vector<pair<mpz_class, int> > factors; //vector<pair<p_i,x_i> >, n = product(p_i^x_i)

    find_previous_factors(n, external_factors); //see if we can find any previously found factors

    for (; n != 1; ++index) {
        bool max_digits_reached = false;
        bool require_both_digits_and_cofactor = both_digits_and_cofactor && max_digits > 0 && max_cofactor > 0;
        int curr_max_cofactor = 0;

        // do trial division here to find small factors, so we can identify the driver
        // this is a bit of a waste since factor() also does it, but it's too messy to try to share the code
        mpz_class n_tmp = n;
        int sm_fac[128];
        memset(sm_fac, 0, sizeof (sm_fac));
        for (size_t j = 0; j < trial_primes.size() && trial_primes[j] <= 127; ++j) {
            while (mpz_divisible_ui_p(n_tmp.get_mpz_t(), trial_primes[j])) {
                sm_fac[ trial_primes[j] ]++;
                mpz_divexact_ui(n_tmp.get_mpz_t(), n_tmp.get_mpz_t(), trial_primes[j]);
            }
        }

        string driver = "No driver";
        bool found_driver;

        // Real drivers will set found_driver to true, and beneficial patterns will set it to false.
        // For guides, we set the default here so that if we're below the digit cutoff we'll keep going,
        // and above it we'll stop.
        if (max_digits > 0 && (int) n.get_str().size() >= max_digits + no_driver_extra)
            found_driver = true;
        else
            found_driver = false;

        if (sm_fac[2] == 1 && sm_fac[3] == 1) {
            driver = "Driver: 2 * 3";
            found_driver = true;
        } else if (sm_fac[2] == 2 && sm_fac[7] == 1) {
            driver = "Driver: 2^2 * 7";
            found_driver = true;
        } else if (sm_fac[2] == 4 && sm_fac[31] == 1) {
            driver = "Driver: 2^4 * 31";
            found_driver = true;
        } else if (sm_fac[2] == 6 && sm_fac[127] == 1) {
            driver = "Driver: 2^6 * 127";
            found_driver = true;
        } else if (sm_fac[2] == 3 && sm_fac[3] == 1 &&
                sm_fac[5] == 1) {
            driver = "Driver: 2^3 * 3 * 5";
            found_driver = true;
        } else if (sm_fac[2] == 3 && sm_fac[3] == 1 &&
                sm_fac[7] == 1) {
            driver = "Driver: 2^3 * 3 * 7";
            found_driver = true;
        } else if (sm_fac[2] == 9 && sm_fac[3] == 1 &&
                sm_fac[11] == 1 && sm_fac[31] == 1) {
            driver = "Driver: 2^9 * 3 * 11 * 31";
            found_driver = true;
        } else if (sm_fac[2] == 3 && sm_fac[3] == 1) {
            driver = "Driver: 2^3 * 3";
            found_driver = true;
        } else if (sm_fac[2] == 2 && sm_fac[3] != 0) {
            driver = "Guide: 2^2 with 3";
        } else if (sm_fac[2] == 3 && sm_fac[5] == 1) {
            driver = "Guide: 2^3 * 5";
        } else if (sm_fac[2] == 3 && sm_fac[7] == 1) {
            driver = "Guide: 2^3 * 7";
        } else if (sm_fac[2] == 3 && sm_fac[3] != 0) {
            driver = "Guide: 2^3 with 3";
        } else if (sm_fac[2] == 5 && sm_fac[3] == 1) {
            driver = "Guide: 2^5 * 3";
        } else if (sm_fac[2] == 1 && sm_fac[3] == 0) {
            driver = "Downdriver!";
            found_driver = false;
        } else if (sm_fac[3] == 0) {
            driver = "No 3";
            found_driver = false;
        } else if (sm_fac[2] >= 8) {
            driver = "Large power of 2";
            found_driver = false;
        }

        // If the sequence decreased, I don't care what pattern we found -- keep going!
        if (n < last_n) {
            driver = "Decreased!";
            found_driver = false;
        }

        cout << "next driver: " << driver << "\r";

        log_msg("\n\n\n", false);
        if (max_digits > 0 && (int) n.get_str().size() >= max_digits) {
            if (require_both_digits_and_cofactor || (!found_driver && exit_only_with_driver)) {
                max_digits_reached = true;
            } else {
                log_msg("*** Max digits " + tostring(max_digits) + " reached for " + seq.get_str() + ":" + tostring(index) +
                        " = " + n.get_str() + " (" + tostring(n.get_str().size()) + " digits)\n");
                cout << endl << "Max digits reached" << endl;
                exit(0);
            }
        }

        log_msg("*** Starting " + seq.get_str() + ":" + tostring(index) + " = " + n.get_str() + " (" + tostring(n.get_str().size()) + " digits)\n");
        curr_max_cofactor = max_cofactor;
        if ((require_both_digits_and_cofactor && !max_digits_reached) || (!found_driver && exit_only_with_driver)) {
            //For this case we don't have both so want a zero max_cofactor to do full factoring
            curr_max_cofactor = 0;
        }
        // factor the current index and save the return value
        bool res_factor = factor(n, factors, external_factors, true, curr_max_cofactor, max_ecm_level);

        merge_factors(factors);

        // save the current index to file even if it is incomplete
        string msg1 = tostring(index) + " .\t ";
        string msg2 = n.get_str() + " = ";
        for (size_t j = 0; j < factors.size(); ++j) {
            if (j) msg2 += " * ";
            msg2 += factors[j].first.get_str() + (factors[j].second > 1 ? ("^" + tostring(factors[j].second)) : "");
        }
        string msg3 = " : " + driver;
        save_result(seq, msg1 + msg2 + "\n");
        cout << msg1 << (factors.size() == 1 && factors[0].second == 1 ? "prp" : "c") << n.get_str().size() << " = " << msg2 << msg3 << endl;

        // look at the return value from factor again we may have to exit
        if (!res_factor) {
            if (require_both_digits_and_cofactor) {
                log_msg("*** Max cofactor size " + tostring(max_cofactor) + " and max digits " + tostring(max_digits) +
                        " reached for sequence " + seq.get_str() + ":" + tostring(index) + " = " + n.get_str()
                        + " (" + tostring(n.get_str().size()) + " digits)\n");
                cout << endl << "Max digits and cofactor reached" << endl;
            } else {
                log_msg("*** Max cofactor size " + tostring(max_cofactor) + " reached for sequence " + seq.get_str()
                        + ":" + tostring(index) + " = " + n.get_str() + " (" + tostring(n.get_str().size()) + " digits)\n");
                cout << endl << "Max cofactor reached" << endl;
            }
            exit(0);
        }

        mpz_class s, verify_n;
        sigma(factors, s, verify_n); //calculate sigma(n) and verify_n = product(factors)

        if (n != verify_n) {
            log_and_print("ERROR: product(factors) != value\n");
            exit(1);
        }

        add_and_check_cycle(index, n);
        check_merge(seq, n, index);

        last_n = n;
        n = s - n; //next sequence value

        skip_ecm = false; //reset flag to do ecm again

        if (quit_after_first) {
            exit(0);
        }
    }

    log_and_print("Woah! Sequence actually ended. It's probably just a bug in the code. ;)\n");

    return 0;
}
