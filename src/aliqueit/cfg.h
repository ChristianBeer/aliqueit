
#ifndef __cfg_h
#define __cfg_h

#include <string>
using namespace std;

class cfg_t {
public:
	bool use_msieve_polyfind;
	bool detect_merge;
    bool stop_failed_gnfs;
	unsigned int trial_cutoff;
	unsigned int gnfs_cutoff;
	bool prefer_yafu;
	bool use_ecmpy;
	string ecm_cmd;
	string ecmpy_cmd;
	string yafu_cmd;
	string msieve_cmd;
	string ggnfs_cmd;
	string ggnfs_clean_cmd;
	unsigned int hide_output_limit;
	unsigned int neat_factor_limit_ecm;
	unsigned int neat_factor_limit_pm1;
	unsigned int neat_factor_limit_pp1;
	string null_device;
	string makedir_cmd;

	string result_file_prefix;
	string log_file;
	bool save_log;
	string ecm_tempfile;

	string ecm_depth;
	unsigned int big_ecm_cutoff;
	float qs_k;
	float qs_m;
	float gnfs_k;
	float gnfs_m;
	float b1scale_ecm;
	float b2scale_ecm;
	float b1scale_pm1;
	float b2scale_pm1;
	float b1scale_pp1;
	float b2scale_pp1;

	int threads;

	cfg_t();
	void read_config_file();
};

extern cfg_t cfg;

#endif //__cfg_h
