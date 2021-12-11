#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hmm.h"
#include <assert.h>
#define MAX_num_obs_seq 10000 // MAX_observation_sequences
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

void forward (HMM *hmm, char *obs_seq, double *alpha, int obs_sqe_len){ // calculate alpha
	int idx = 0;
    // init
	for (int n = 0 ; n < hmm->state_num; n++){
        idx = n; // 0*(hmm->state_num)+n
		alpha[idx] = hmm->initial[n] * hmm->observation[obs_seq[0] - 'A'][n];
        // printf("**alpha[%d] = %f\n", idx, alpha[idx]);
	}

    double sum = 0.;
	// induction 
    for (int t = 1; t < obs_sqe_len; t++) { 
        for (int j = 0; j < hmm->state_num; j++) {  // for all states
            sum = 0.;
            for (int i = 0; i < hmm->state_num; i++) { // accumulate the states from the previous T 
                idx = (t-1)*(hmm->state_num) + i;
                sum += ((alpha[idx])*(hmm->transition[i][j]));//part; // sum += (alpha[abgIDX((t-1), i)])*(hmm->transition[i][j]); 會死去不知道為啥QQ
            }
            idx = t*(hmm->state_num) + j;
            alpha[idx] = sum * hmm->observation[obs_seq[t] - 'A'][j];
            // printf("alpha[%d][%d] = %f\n", t, j, alpha[idx]);
        }
    } 
}

void backward (HMM *hmm, char *obs_seq, double *beta, int obs_sqe_len){ // calculate beta
	// init
    int idx = 0;
    for (int n = 0; n < hmm->state_num; n++) {
        // beta[obs_sqe_len - 1][n] = 1; 
        idx = ((obs_sqe_len - 1)*(hmm->state_num) + (n));
        beta[idx] = 1;
        // printf("idx = %d, beta[idx] = %f\n", idx, beta[idx] );
    }

    double sum = 0.;
	// induction 
    for (int t = obs_sqe_len-2; t >= 0; t--) {
        for (int i = 0; i < hmm->state_num; i++) { // for the previous i states @ T = (t-1)
            sum = 0.;
            for (int j = 0; j < hmm->state_num; j++) { // from the j states to the i state @ T = (t)
                idx = (t+1)*(hmm->state_num) + j;
                // sum += hmm->transition[i][j] * hmm->observation[obs_seq[t+1] - 'A'][j] * beta[t+1][j];
                sum += ((hmm->transition[i][j])*(hmm->observation[obs_seq[t+1] - 'A'][j])*(beta[idx]));
            }												   
            // beta[t][i] = sum;
            idx = t*(hmm->state_num) + i;
            beta[idx] = sum;
            // printf("beta[%d][%d] = %f\n", t, i, beta[idx]);
        } 
    }

}

void calculate_gamma_epsilon (HMM *hmm, char *obs_seq, double *alpha, double *beta, double *gamma, double *epsilon, double *prob_ob, int obs_sqe_len, double *gamma_mom) {
	// P(O_bar, qt = i|lambda):  prob_ob += alpha[][]*beta[][];
    int idx_a = 0;
    // gamma的方母
    for (int t = 0; t < obs_sqe_len; t++) {
		prob_ob[t] = 0;
        for (int n = 0; n < (hmm->state_num); n++){
            idx_a = t*(hmm->state_num) + n;
			prob_ob[t] += alpha[idx_a]*beta[idx_a];
		}
	}

    int idx_b = 0, idx_e = 0;
	// // // gamma, 分母
 //    for (int t = 0; t < obs_sqe_len; t++) {
 //        for (int n = 0; n < hmm->state_num; n++) {
 //            idx_a = t*(hmm->state_num) + n;
 //            gamma_mom[idx_a] = 0;
 //            for (int j = 0; j < (hmm->state_num); j++){
 //                idx_b = t*(hmm->state_num) + j;
 //                gamma_mom[idx_a] += (alpha[idx_b] * beta[idx_b]);
 //            }
 //        // printf("prob_ob[%d] = %f,  gamma_mom[%d] = %f\n", t, prob_ob[t],  idx_a, gamma_mom[idx_a]);
 //        }
 //    }


    // gamma
    for (int t = 0; t < obs_sqe_len; t++) {
        for (int n = 0; n < hmm->state_num; n++) {
            idx_a = t*(hmm->state_num) + n;
            // gamma ＝ 分子/分母
            // gamma[idx_a] = (alpha[idx_a] * beta[idx_a]) / gamma_mom[idx_a];
            gamma[idx_a] = (alpha[idx_a] * beta[idx_a]) / prob_ob[t];
        }
    }
    // epsilon
    // double epsilon_mom[(MAX_obs_sqe_len-1)][MAX_STATE][MAX_STATE] = {0.};
    for (int t = 0; t < obs_sqe_len-1; t++) {
        for (int i = 0; i < hmm->state_num; i++) {
            // 分子
            for (int j = 0; j < hmm->state_num; j++) {
                idx_a = t*(hmm->state_num) + i;
                idx_b = (t+1)*(hmm->state_num) + j;
                idx_e = t*(hmm->state_num)*(hmm->state_num) + i*(hmm->state_num) + j;
                // epsilon[idx_e] = ((alpha[idx_a]) * (hmm->transition[i][j]) * (hmm->observation[obs_seq[t+1] - 'A'][j]) * (beta[idx_b])) / (gamma_mom[idx_a]);
                epsilon[idx_e] = ((alpha[idx_a]) * (hmm->transition[i][j]) * (hmm->observation[obs_seq[t+1] - 'A'][j]) * (beta[idx_b])) / prob_ob[t];
            }
        }
    }
}

//Expected: alpha, beta, gamma, epsilon
void E_step (HMM *hmm, char *obs_seq, double *alpha, double *beta, double *gamma, double *epsilon, int num_obs_seq, int obs_sqe_len, double *prob_ob, double *gamma_mom) { 
	forward (hmm, obs_seq, alpha, obs_sqe_len);
	backward (hmm, obs_seq, beta, obs_sqe_len);
	calculate_gamma_epsilon (hmm, obs_seq, alpha, beta, gamma, epsilon, prob_ob, obs_sqe_len, gamma_mom);	
}

void M_step (HMM *hmm, char *obs_seq, double *alpha, double *beta, double *gamma, double *epsilon, double *prob_ob, int obs_sqe_len, double p[MAX_STATE][2], double a[MAX_STATE][MAX_STATE][2], double b[MAX_STATE][MAX_STATE][2]) {
    double sum, part;
    int idx = 0, idx_2 = 0;
    // accumulate pi
    for (int i = 0; i < (hmm->state_num); i++) {
        idx = i; // gamma idx = 0*MAX_STATE + i
        p[i][0] += gamma[idx];
    }

    // accumulate A
    idx = 0;
    // 分子
    for (int i = 0; i < (hmm->state_num); i++) {
        for (int j = 0; j < (hmm->state_num); j++) {
            part = 0;
            for (int t = 0; t < obs_sqe_len-1; t++) {
                idx = t*(hmm->state_num)*(hmm->state_num) + i*(hmm->state_num) + j;
                part += epsilon[idx];
            }
            a[i][j][0] += part; 
        }
    }
    // 分母
    idx = 0;
    for (int i = 0; i < (hmm->state_num); i++) {
        for (int j = 0; j < (hmm->state_num); j++) {
        	sum = 0;
            for (int t = 0; t < obs_sqe_len-1; t++) {
                idx = t*(hmm->state_num) + i;
                sum += gamma[idx];
            }
        	a[i][j][1] += sum; // 分母
        }
    }

    // accumulate B
    for (int i = 0; i < (hmm->state_num); i++) {
        for (int s = 0; s < (hmm->state_num); s++) {
            sum = 0, part = 0;
            for (int t = 0; t < obs_sqe_len; t++) {
                idx = t*(hmm->state_num) + i;
                if ((obs_seq[t] - 'A') == s) {
                    part += gamma[idx];                
                }
                sum += gamma[idx];
            }
            b[s][i][0] += part;
            b[s][i][1] += sum;
        }
    }
}

// Re-estimate: "update" alpha, beta with gamma and epsilon 
void update (HMM *hmm , double p[MAX_STATE][2], double a[MAX_STATE][MAX_STATE][2], double b[MAX_STATE][MAX_STATE][2], int num_obs_seq) {
    // update pi, A, B
	for (int n = 0 ; n < hmm->state_num ; n++){
        hmm->initial[n] = p[n][0]/num_obs_seq;
        // hmm->initial[n] = p[n][0]/p[n][1];
        // printf("p[%d][1] = %f, num_obs_seq = %d\n", n, p[n][1], num_obs_seq);
	}

    for (int i = 0; i < (hmm->state_num); i++) {
        for (int j = 0; j < (hmm->state_num); j++) {
            hmm->transition[i][j]  = a[i][j][0] / a[i][j][1];
            // if (b[i][j][0] / b[i][j][1] > 0){
            //     printf("b[%d][%d][0]  = %f,  b[%d][%d][1] = %f\n", i, j, b[i][j][0], i, j, b[i][j][1]);
            // }
            hmm->observation[i][j] = b[i][j][0] / b[i][j][1];
        }
    }
}

void train (HMM *hmm , int iteration, char *obs_seq, int num_obs_seq, int obs_sqe_len, double *alpha, double *beta, double *gamma, double *epsilon, double *prob_ob, double *gamma_mom){
// Baum-Welch algorithm: A generalized expectation-maximization (EM) algorithm
	for (int i = 0; i < iteration; i++) {
        double p[MAX_obs_sqe_len][2] = {0};
        double a[MAX_STATE][MAX_STATE][2] = {0}; // 分子、分母
        double b[MAX_STATE][MAX_STATE][2] = {0}; // 分子、分母
		for (int j = 0; j < num_obs_seq; j++) {
            // printf("obs_seq = %s\n", obs_seq + j*obs_sqe_len);
			E_step (hmm, obs_seq + j*obs_sqe_len, alpha, beta, gamma, epsilon, num_obs_seq, obs_sqe_len, prob_ob, gamma_mom);//obs_seq start idx: j*obs_sqe_len
			M_step (hmm, obs_seq + j*obs_sqe_len, alpha, beta, gamma, epsilon, prob_ob, obs_sqe_len, p, a, b);
		}
		update (hmm, p, a, b, num_obs_seq);
	}
}

int main(int argc, char **argv) {
    if (argc < 5) {
        printf("usage example:\n"); 
        printf("./train iteration model_init.txt ./data/train_seq_01.txt model_01.txt\n");
        exit(0);
    }

    // argv: (1, iteration), (2, model_init.txt), (3, model_seq), (4, model_output)

	HMM hmm;
	loadHMM (&hmm, argv[2]);  // initial model
	
	char *obs_seq = (char*)malloc(MAX_num_obs_seq*MAX_obs_sqe_len*sizeof(char));
	memset (obs_seq, 0, MAX_num_obs_seq*MAX_obs_sqe_len*sizeof(char));
	
	double *alpha = (double*)malloc(MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	double *beta  = (double*)malloc(MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	double *gamma = (double*)malloc(MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	double *epsilon = (double*)malloc((MAX_obs_sqe_len-1)*MAX_STATE*MAX_STATE*sizeof(double));
	memset (alpha, 0, MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	memset (beta, 0, MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	memset (gamma, 0, MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	memset (epsilon, 0, (MAX_obs_sqe_len-1)*MAX_STATE*MAX_STATE*sizeof(double));
	
	double *prob_ob = (double*)malloc(MAX_obs_sqe_len*sizeof(double)); // P(O_bar|lambda), update pi needed
	memset (prob_ob, 0, MAX_obs_sqe_len*sizeof(double));
    double *gamma_mom = (double*)malloc(MAX_obs_sqe_len*MAX_STATE*sizeof(double));
    memset (gamma_mom, 0, MAX_obs_sqe_len*MAX_STATE*sizeof(double));
	
    int num_obs_seq = 0, obs_sqe_len = 0;
	int iteration = atoi(argv[1]);

	read_file (argv[3], obs_seq, &num_obs_seq, &obs_sqe_len);
	// printf("after: num_obs_seq = %d, obs_sqe_len = %d\n", num_obs_seq, obs_sqe_len);

	train (&hmm, iteration, obs_seq, num_obs_seq, obs_sqe_len, alpha, beta, gamma, epsilon, prob_ob, gamma_mom);
	

	FILE *Output = fopen(argv[4] , "w");
	dumpHMM(Output, &hmm); // Output
    fclose(Output);
    return 0;
}

