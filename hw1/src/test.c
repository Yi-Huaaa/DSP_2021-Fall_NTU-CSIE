#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"
#include <assert.h>
#define num_model 5 // max: 5
#define MAX_num_obs_seq 2500
#define MAX_obs_sqe_len 51

void print_obs_seq (char *obs_seq, int num_obs_seq, int obs_sqe_len) {
	printf("print_obs_seq:\n");
	for (int i = 0; i < num_obs_seq; i++) {
		for (int j = 0; j < obs_sqe_len; j++) {
			printf("%c", obs_seq[i*obs_sqe_len + j]);
		}
		printf("\n");
	}
	printf("---------------\n");
}

void read_file (char *file_name, char *obs_seq, int *num_obs_seq, int *obs_sqe_len) {	
	FILE *instance = fopen(file_name , "r");
    assert(instance != NULL);

    fscanf(instance , "%s" , &obs_seq[0]);
    (*obs_sqe_len) = strlen(obs_seq);
    (*num_obs_seq) ++;

    while (!feof(instance)) {
		fscanf(instance , "%s" , &obs_seq[(*num_obs_seq)*(*obs_sqe_len)]);
		(*num_obs_seq) ++;
    }
    fclose(instance);
    (*num_obs_seq) --;

	// print_obs_seq (obs_seq, (*num_obs_seq), (*obs_sqe_len));
	// printf("num_obs_seq = %d, obs_sqe_len = %d\n", (*num_obs_seq), (*obs_sqe_len));
}

double Viterbi (HMM *hmm, char *obs_seq, int obs_sqe_len, double *delta){ // m: model number
	int idx = 0;
	double tmp_MAX = 0., result = 0.;
	//init
	for (int n = 0; n < hmm->state_num; n++) {
		idx = n; // 0*(hmm->state_num) + n
		delta[idx] = hmm->initial[n] * hmm->observation[obs_seq[0] - 'A'][n];
		// printf("delta[%d] = %f\n", idx, delta[idx]);

	}
	// Recursion
	for (int t = 1; t < obs_sqe_len; t++) {
		for (int j = 0; j < hmm->state_num; j++) {
			tmp_MAX = 0.;
			for (int i = 0; i < hmm->state_num; i++) {
				idx = (t-1)*(hmm->state_num) + i;
				result = delta[idx] * hmm->transition[i][j];
				// printf("delta[%d] = %f, hmm->transition[%d][%d] = %f, result = %f\n", idx, delta[idx], i, j, hmm->transition[i][j], result);
				if (tmp_MAX < result) {
					tmp_MAX = result;
				}
			}
			tmp_MAX *= hmm->observation[obs_seq[t] - 'A'][j];
			idx = (t)*(hmm->state_num) + j;
			delta[idx] = tmp_MAX;
			// printf("delta[%d] = %f\n", idx, delta[idx]);

		}
	}
	tmp_MAX = 0;
    for (int i = 0; i < (hmm->state_num); i++) {
    	idx = (obs_sqe_len-1)*(hmm->state_num) + i; // [obs_sqe_len-1][i]
        if (tmp_MAX < delta[idx]) {
        	tmp_MAX = delta[idx];
        }
    }
    // printf("tmp_MAX = %f\n", tmp_MAX);
	return tmp_MAX;
}

void test (HMM hmm[num_model], char *obs_seq_file, char *result_file, char *obs_seq, int num_obs_seq, int obs_sqe_len, double *delta){
	double max = 0., res = 0.;
    int who = 0;
	read_file (obs_seq_file, obs_seq, &num_obs_seq, &obs_sqe_len);
	
	FILE *dump_file = fopen(result_file , "w");

	for (int j = 0 ; j < num_obs_seq ; j++) {
        max = 0., who = 0;		
		for (int m = 0; m < num_model; m++){
			res = Viterbi(&hmm[m], obs_seq + j*obs_sqe_len, obs_sqe_len, delta);
			// printf("res = %f\n", res);
			if (max < res){
                max = res; 
                who = m;
			}
		} 
		// fprintf(dump_file , "%s%d%s %e\n" , "model_0", (who+1), ".txt num" , max);
		fprintf(dump_file , "%s%d%s %e\n" , "model_0", (who+1), ".txt" , max);
	}
	fclose(dump_file);
}


int main(int argc, char **argv) {
    if (argc < 4) {
    	printf("usage example:\n");
    	// printf("./test modellist.txt test_seq.txt result.txt test_lbl.txt");
    	printf("./test modellist.txt test_seq.txt result.txt");
        exit(0);
    }
    // argv: (1, modellist), (2, test_seq.txt), (3, result.txt), (4, test_lbl.txt)

	HMM hmm[num_model];
	load_models(argv[1] , hmm, num_model);

	char *obs_seq = (char*)malloc(MAX_num_obs_seq*MAX_obs_sqe_len*sizeof(char));
	memset (obs_seq, 0, MAX_num_obs_seq*MAX_obs_sqe_len*sizeof(char));

	double *delta = (double*)malloc(MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	memset (delta, 0, MAX_obs_sqe_len*MAX_STATE*sizeof(double));


	int num_obs_seq = 0, obs_sqe_len = 0;
	test (hmm , argv[2] , argv[3], obs_seq, num_obs_seq, obs_sqe_len, delta);
	
// calculate accuracy
 //    char acc_file [8] = "acc.txt";
	// char str1 [500] = {0};
	// char str2 [500] = {0};

 //    if (strcmp(argv[2] , "./data/test_seq.txt") == 0) { 
	// 	FILE *res = fopen(argv[3] , "r");
	// 	FILE *answer = fopen(argv[4] , "r");
	// 	FILE *acc = fopen(acc_file , "w");

	// 	// double max = 0.;
	// 	double max;
	// 	double trueAns = 0, allAnsewer = 0;

	//     while (!feof(res)) {
	// 		fscanf(res, "%s %lf" , &str1, &max);//這兩行有bug
	// 		// printf("str1 = %s\n", str1);
	// 		fscanf(answer, "%s  %lf" , &str2, &max);//這兩行有bug
 
	// 		if (strcmp(str1 , str2) == 0){
	// 			trueAns++;
	// 		}
	// 		allAnsewer++;
	//     }

	// 	fprintf(acc, "%f\n" , trueAns/allAnsewer);

	// 	fclose(res);
	// 	fclose(answer);
	// 	fclose(acc);
 //    }
	return 0;
}
