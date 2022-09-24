
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "cfg.h"
#include "misc.h"
using namespace std;

cfg_t::cfg_t() {
    //------ default config, overridden by config file -------
    //Use ecm.py (for easy ecm multithreading) instead of regular ecm?
    use_ecmpy = false;

    //Use msieve to find a gnfs poly? If false, ggnfs will have to find its own.
    use_msieve_polyfind = true;

    //Stop if sequence merges with an earlier sequence?
    detect_merge = true;

    //Stop on qs or gnfs failures, otherwise try again or run autoincreasing ecm until factor found
    stop_on_failure = false;

    //When -t is used exit with 0 if sequence in elf is correct and terminates, exit nonzero otherwise
    verify_terminations = false;

    //Trial factor with primes < this value.
    trial_cutoff = 10000;

    //Use GNFS for composites >= this many digits.
    gnfs_cutoff = 99;

    //Use yafu before msieve?
    //bool prefer_yafu = false;
    prefer_yafu = true;

    //Name (and path, if needed) of your gmp-ecm executable.
    ecm_cmd = "ecm";

    //Name (and path, if needed) of your ecm.py executable.
    ecmpy_cmd = "python.exe ecm.py";

    //Name (and path, if needed) of your yafu executable.
    yafu_cmd = "yafu";

    //Name (and path, if needed) of your msieve executable.
    msieve_cmd = "msieve";

    //Name (and path, if needed) of your ggnfs executable.
    ggnfs_cmd = "c:\\strawberry\\perl\\bin\\perl.exe \"c:\\program files (x86)\\ggnfs\\bin\\factMsieve.pl\"";

    //Take out the trash if factorisation was successful.
    ggnfs_clean_cmd = "del cols*.* deps*.* factor.easy lpindex*.* rels*.* "
            " numFF.txt numrels.txt tmp.job spmat sp-index test.fb matsave matsave.bak "
            " test.dat*";

    //Don't show msieve/yafu output on screen for composites smaller than this limit (avoids screen spam).
    hide_output_limit = 70;

    //Print a cheerful msg if an ecm factor with at least this many digits is found.
    neat_factor_limit_ecm = 45;

    //Print a cheerful msg if a P-1 factor with at least this many digits is found.
    neat_factor_limit_pm1 = 45;

    //Print a cheerful msg if a P+1 factor with at least this many digits is found.
    neat_factor_limit_pp1 = 40;

    //Null device (for redirecting unwanted output to).
    null_device = "nul"; //windows
    //null_device = "/dev/null";	//linux

    //Shell command to create a directory
    makedir_cmd = "mkdir";



    //Save result output to this file + "<sequence_number>.elf"
    result_file_prefix = "alq_";

    //Save all log output from external programs to this file...
    log_file = "aliqueit.log";
    //  ...if this is true, otherwise discard it
    save_log = true;

    //Temp file for ecm output. Moved to <log_file> after each run.
    ecm_tempfile = "aliqueit_ecm_temp.log";


    //ECM depth: comma-separated number pairs "<digits> <curves>". Used to select the number of curves to run.
    //For a composite of x digits the highest pair whose <digits> element <= x is chosen.
    //[gmp-]ecm is then run with #curves=<curves> ("-c <curves>") and B1 auto-increasing from 1 ("-I 1").
    //20 curves ~ B1=7e3
    //50 curves ~ B1=26e3
    //100 curves ~ B1=78e3
    //200 curves ~ B1=0.23e6
    //300 curves ~ B1=0.44e6
    //400 curves ~ B1=0.65e6
    //500 curves ~ B1=0.88e6
    //600 curves ~ B1=1.12e6
    //700 curves ~ B1=1.38e6
    //800 curves ~ B1=1.65e6
    //1000 curves ~ B1=2.25e6
    //For <65 digits I've experimentally found these values to be good on my machine. Improvements are always welcome!
    //NOTE: These values are only used for composites < <big_ecm_cutoff> digits. See that setting below.
    ecm_depth =
            "45 10, 50 15, 55 20, 60 30, 65 40, 70 60, 73 90, 76 120, 79 150, 82 180, "
            "83 200, 86 300, 90 400, 93 500, 96 600, 99 700, 102 800, 106 900, 110 1000, 114 1100, 118 1200, 122 1300, ";

    //If composite has more than <big_ecm_cutoff> digits, we will try to do ecm for about 1/4 of the estimated time qs/nfs would take.
    //I've used logs of previous factorisations on my machine to construct the estimates that are used. They have been converted to
    //linear formulas taking the composite size as input and giving maximum factor size (see Table 1 in gmp-ecm readme) as output.
    //See also qs_k, qs_m, gnfs_k, and gnfs_m below.
    big_ecm_cutoff = 65;

    //Formulae used to determine the maximum factor size we will do ecm to.
    //For QS: <qs_k> * input_digits + <qs_m>
    qs_k = 0.448f;
    qs_m = -11.26f;

    //For GNFS: <gnfs_k> * input_digits + <gnfs_m>
    gnfs_k = 0.235f;
    gnfs_m = 9.4f;

    //Multiply default B1 values by these constants.
    b1scale_ecm = 1.0f;
    b1scale_pm1 = 1.0f;
    b1scale_pp1 = 1.0f;

    //max number of threads to use
    threads = 1;
    //------ config end -------
}

bool get_bool(string val) {
    val = tolower(val);
    if (val != "false" && val != "true") {
        cout << "ERROR: expected boolean value (true/false), but got: " << val << endl;
        exit(1);
    }
    return val == "true";
}

unsigned int get_uint(string val) {
    if (!isnumber(val)) {
        cout << "ERROR: expected numeric value, but got: " << val << endl;
        exit(1);
    }
    return atoi(val.c_str());
}

float get_float(string val) {
    float f = (float) atof(val.c_str());
    if (val.size() && val[0] == '+' || val[0] == '-') {
        val = val.substr(1); //ok, skip by sign
    }
    if (isnumber(val)) {
        return f;
    }
    size_t o = val.find('.');
    if (o == val.npos || !isnumber(val.substr(0, o)) || !isnumber(val.substr(o + 1))) {
        cout << "ERROR: expected numeric value, but got: " << val << endl;
        exit(1);
    }
    return f;
}

void cfg_t::read_config_file() {
    ifstream f("aliqueit.ini");
    if (!f.is_open()) {
        cout << "WARNING: couldn't open config file." << endl;
        return;
    }
    cout << "Reading config file..." << endl;
    string line;
    while (getline(f, line)) {
        line = trim(line);
        if (!line.size() || (line[0] == '#' || line[0] == ';' || line[0] == '/')) continue; //comment
        if (line[0] == '[') continue; //section name. unused...
        size_t o;
        if ((o = line.find("=")) == line.npos) {
            cout << "ERROR: unrecognized line: " << line << endl;
            exit(1);
        }
        string arg = tolower(trim(line.substr(0, o))); //arg name
        string val = trim(line.substr(o + 1)); //value
        //cout << arg << " = " << val << endl;

        //assign values...
        if (arg == "trial_cutoff") {
            trial_cutoff = get_uint(val);
        } else if (arg == "gnfs_cutoff") {
            gnfs_cutoff = get_uint(val);
        } else if (arg == "prefer_yafu") {
            prefer_yafu = get_bool(val);
        } else if (arg == "ecm_tempfile") {
            ecm_tempfile = val;
        } else if (arg == "ecm_cmd") {
            ecm_cmd = val;
        } else if (arg == "ecmpy_cmd") {
            ecmpy_cmd = val;
        } else if (arg == "yafu_cmd") {
            yafu_cmd = val;
        } else if (arg == "msieve_cmd") {
            msieve_cmd = val;
        } else if (arg == "ggnfs_cmd") {
            ggnfs_cmd = val;
        } else if (arg == "ggnfs_clean_cmd") {
            ggnfs_clean_cmd = val;
        } else if (arg == "result_file_prefix") {
            result_file_prefix = val;
        } else if (arg == "log_file") {
            log_file = val;
        } else if (arg == "save_log") {
            save_log = get_bool(val);
        } else if (arg == "makedir_cmd") {
            makedir_cmd = val;
        } else if (arg == "ecm_depth") {
            ecm_depth = val;
        } else if (arg == "hide_output_limit") {
            hide_output_limit = get_uint(val);
        } else if (arg == "null_device") {
            null_device = val;
        } else if (arg == "detect_merge") {
            detect_merge = get_bool(val);
        } else if (arg == "stop_on_failure") {
            stop_on_failure = get_bool(val);
        } else if (arg == "verify_terminations") {
            verify_terminations = get_bool(val);
        } else if (arg == "big_ecm_cutoff") {
            big_ecm_cutoff = get_uint(val);
        } else if (arg == "neat_factor_limit_ecm") {
            neat_factor_limit_ecm = get_uint(val);
        } else if (arg == "neat_factor_limit_pm1") {
            neat_factor_limit_pm1 = get_uint(val);
        } else if (arg == "neat_factor_limit_pp1") {
            neat_factor_limit_pp1 = get_uint(val);
        } else if (arg == "qs_k") {
            qs_k = get_float(val);
        } else if (arg == "qs_m") {
            qs_m = get_float(val);
        } else if (arg == "gnfs_k") {
            gnfs_k = get_float(val);
        } else if (arg == "gnfs_m") {
            gnfs_m = get_float(val);
        } else if (arg == "use_msieve_polyfind") {
            use_msieve_polyfind = get_bool(val);
        } else if (arg == "use_ecmpy") {
            use_ecmpy = get_bool(val);
        } else if (arg == "b1scale_ecm") {
            b1scale_ecm = get_float(val);
        } else if (arg == "b1scale_pm1") {
            b1scale_pm1 = get_float(val);
        } else if (arg == "b1scale_pp1") {
            b1scale_pp1 = get_float(val);
        } else if (arg == "threads") {
            threads = get_uint(val);
        } else {
            cout << "ERROR: unrecognized arg: " << arg << endl;
            exit(1);
        }
    }
    f.close();
}


cfg_t cfg;
