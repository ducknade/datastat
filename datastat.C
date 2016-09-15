#include <stdio.h>
#include <readline/readline.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <string>

#include "statjks.h"

using namespace std;

FILE *fin = stdin;
const char *delim = " \t";
bool show_avg = true;
bool show_dev = false;
bool show_min = false;
bool show_max = false;
bool show_cnt = false;
bool show_sum = false;
bool show_atc = false;
bool show_jkn = false;

struct record {
	vector<double> v_sum;
	vector<double> v_sqr;
	vector<double> v_min;
	vector<double> v_max;
	vector<long> v_num;
	vector< vector<double> > v_val;
	record(): v_sum(), v_sqr(), v_min(), v_max(), v_num(), v_val() {}
};

void usage() {
	printf("Usage: datastat [options] [filename]\n");
	printf("  Options:\n");
	printf("    -h|--help ....... This help message\n");
	printf("    --no-avg ........ Suppress average\n");
	printf("    --dev ........... Show standard deviation\n");
	printf("    --min ........... Show minimum\n");
	printf("    --max ........... Show maximum\n");
	printf("    --sum ........... Show sum\n");
	printf("    --cnt ........... Show count of values\n");
	printf("    --atc ........... Show autocorrelation time in detail\n");
	printf("    --jkn ........... Show jackknife standard deviation\n");
}

int ceil_tol(double input, double tol){
    int ceil = (int)std::ceil(input);
    int trun = (int)input;
    if(input - trun < 1e-4){
        return trun;
    }else{
        return ceil;
    }    
}

int main(int argc, char *argv[]){
	--argc;  
	++argv;
	while(argc > 0){
		if(strcmp(*argv, "-h") == 0 || strcmp(*argv, "--help") == 0){
			usage();
			exit(0);
		}else if(strcmp(*argv, "--no-avg") == 0){
		 	show_avg = false;
		}else if(strcmp(*argv, "--dev") == 0){
		 	show_dev = true;
		}else if(strcmp(*argv, "--min") == 0){
		 	show_min = true;
		}else if(strcmp(*argv, "--max") == 0){
		 	show_max = true;
		}else if(strcmp(*argv, "--sum") == 0){
		 	show_sum = true;
		}else if(strcmp(*argv, "--cnt") == 0){
		 	show_cnt = true;
		}else if(strcmp(*argv, "--atc") == 0){ // auto correlation
		 	show_atc = true;
		}else if(strcmp(*argv, "--jkn") == 0){ // jackknife 
		 	show_jkn = true;
		}else{
		  	fin = fopen(*argv, "r");
		}
		--argc;
		++argv;
	}

	record accum;
	
	while (!feof(fin)) {
		char *line = NULL;
		size_t line_size = 0;
		ssize_t rv;
		rv = getline(&line, &line_size, fin);
		if (rv < 0 || line == NULL)
			break;
		if (line[0] == '#')	// Can have comments in the file
		  	continue;
		if (line[strlen(line) - 1] == '\n')
		    line[strlen(line) - 1] = '\0';
		
		vector<string> values;	// vector of token values for current line
		char *tok = strtok(line, delim);
		while(tok != NULL){
			values.push_back(string(tok));
		  	tok = strtok(NULL, delim);
		}
		free(line);

		for (int i = 0; i < (int)values.size(); ++i) {
			const char *s = values[i].c_str();
			double d;
			sscanf(s, "%lf", &d);                         
				// d now contains the double value of the read string
			if (i >= (int)accum.v_sum.size()) {
				accum.v_sum.push_back(d);
				accum.v_sqr.push_back(d*d);
				accum.v_min.push_back(d);
				accum.v_max.push_back(d);
				accum.v_val.push_back(vector<double>()); 
				accum.v_val[i].push_back(d);
				accum.v_num.push_back(1);
			} else {
				accum.v_sum[i] += d;
				accum.v_sqr[i] += d*d;
				accum.v_min[i] = min(accum.v_min[i], d);
				accum.v_max[i] = max(accum.v_max[i], d);
				accum.v_val[i].push_back(d);   
				accum.v_num[i]++;
			}
		}
	}
	// end of reading and processing file... now output final info
	printf("#");
	if (show_avg) {
	   printf("avg\t\t");
	}
	if (show_dev) {
	   printf("dev\t\t");
	}
	if (show_min) {
	   printf("min\t\t");
	}
	if (show_max) {
	   printf("max\t\t");
	}
	if (show_sum) {
	   printf("sum\t\t");
	}
	if(show_jkn){
	   printf("atc\t\tjkn\t\t");
	}
	if(show_cnt){
	   printf("cnt");
	}
	printf("\n");

	for (int i = 0; i < (int)accum.v_sum.size(); ++i) {
		double avg = accum.v_sum[i] / accum.v_val[i].size();
		if (show_avg) {
			printf("%.6e\t", avg);
		}
		if (show_dev) {
			printf("%.6e\t", 
			sqrt(accum.v_sqr[i] / accum.v_val[i].size() - avg * avg));
		}
		if (show_min) {
			printf("%.6e\t", accum.v_min[i]);
		}
		if (show_max) {
			printf("%.6e\t", accum.v_max[i]);
		}
		if (show_sum) {
			printf("%.6e\t", accum.v_sum[i]);
		}

		// start jackson
		if(show_jkn){
			double atc = autoCorrelation(accum.v_val[i]);
			int binSize = ceil_tol(atc);
			double avg_;
			double jkn = jackknife(accum.v_val[i].data(),
					accum.v_val[i].size(), binSize, avg_);
			printf("%.6e\t", atc);
			printf("%.6e\t", jkn);
		}
		// end jackson
		
		if (show_cnt) {
		printf("%lu", accum.v_num[i]);
		}

		if(show_atc){
			vector<double> series;
			double atc = autoCorrelation(accum.v_val[i], series);
			printf("\n---- auto-correlation time accumulate ----\n");
			for(int j = 0; j < series.size(); j++){
				printf("%d\t%.3f\n", j, series[j]);
			}
			printf("atc\t%.3f\n", atc);
			printf("---- auto-correlation time accumulate ----\n");
		}

		printf("\n");
	}
}
