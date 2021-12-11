#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"
#include <assert.h>
#define num_model 5 // max: 5
#define MAX_num_obs_seq 2500
#define MAX_obs_sqe_len 51

int main(int argc, char **argv) {
    if (argc < 4) {
    	printf("usage example:\n");
    	printf("./test modellist.txt test_seq.txt result.txt test_lbl.txt");
        exit(0);
    }
    // argv: (1, modellist), (2, test_seq.txt), (3, result.txt), (4, test_lbl.txt)

    char acc_file [8] = "acc.txt";
	char str1 [500] = {0};
	char str2 [500] = {0};

    if (strcmp(argv[2] , "./data/test_seq.txt") == 0) { 
		FILE *res = fopen(argv[3] , "r");
		FILE *answer = fopen(argv[4] , "r");
		FILE *acc = fopen(acc_file , "w");

		// double max = 0.;
		double max;
		double trueAns = 0, allAnsewer = 0;

	    while (!feof(res)) {
			fscanf(res, "%s %lf" , &str1, &max);
			// printf("str1 = %s\n", str1);
			fscanf(answer, "%s  %lf" , &str2, &max);
 
			if (strcmp(str1 , str2) == 0){
				trueAns++;
			}
			allAnsewer++;
	    }

		fprintf(acc, "%f\n" , trueAns/allAnsewer);

		fclose(res);
		fclose(answer);
		fclose(acc);
    }
	return 0;
}