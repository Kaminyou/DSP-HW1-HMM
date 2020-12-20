#include "hmm.h"
#include <iostream>
#include <fstream>
#include <string.h>
using namespace std;

#define TOTAL_STATE 6
#define TOTAL_OBV 6
#define SEQ_LEN 100

typedef struct{
	double alpha[SEQ_LEN][TOTAL_STATE];
	double beta[SEQ_LEN][TOTAL_STATE];
	double gamma[SEQ_LEN][TOTAL_STATE];
	double epsilon[SEQ_LEN][TOTAL_STATE][TOTAL_STATE];
} Temp_Paras;

typedef struct{
	double pi[TOTAL_STATE];
	double a_u[TOTAL_STATE][TOTAL_STATE];
	double a_l[TOTAL_STATE];
	double b_u[TOTAL_STATE][TOTAL_OBV];
    double b_l[TOTAL_STATE];
    int times;
} Sum_Paras;

void calc_paras(HMM *, Temp_Paras *, Sum_Paras *, string sequence);
void update_hmm(HMM *, Sum_Paras *);

int main(int argc, char **argv){
    
	int iter = stoi(argv[1]);
	const char* init_hmm = argv[2];
	const char* input_file = argv[3];
	const char* output_file = argv[4];
	
	HMM hmm;
    loadHMM(&hmm, init_hmm);

	Temp_Paras temp_paras;
    Sum_Paras sum_paras;
    string sequence="";
    
	for (int i = 0; i < iter; i++){

        for (int n_i=0; n_i<TOTAL_STATE; n_i++){
            for (int n_j=0; n_j<TOTAL_STATE; n_j++){
                sum_paras.a_u[n_i][n_j] = 0.0;
            }
            for (int k=0; k<TOTAL_OBV; k++){
                sum_paras.b_u[n_i][k] = 0.0;
            }
            sum_paras.a_l[n_i] = 0.0;
            sum_paras.b_l[n_i] = 0.0;
            sum_paras.pi[n_i] = 0.0;
        }
        sum_paras.times = 0;


		ifstream ifs(input_file, ifstream::in);
		while(getline(ifs, sequence)){
            calc_paras(&hmm, &temp_paras, &sum_paras, sequence);
		}
		update_hmm(&hmm, &sum_paras);
		ifs.close();
	}
	dumpHMM(open_or_die(output_file, "w"), &hmm);
	return 0;
}


void calc_paras(HMM *hmm, Temp_Paras *temp_paras, Sum_Paras *sum_paras, string sequence){
    int sequence_len = sequence.size();

	// alpha t: T n: N
    for (int n = 0; n < TOTAL_STATE; n++){
		temp_paras->alpha[0][n] = (hmm->initial[n]) * (hmm->observation[sequence[0]-65][n]);
	}
    double pre_alpha_sum;
	for (int t = 1; t < sequence_len; t++){
		for (int n = 0; n < TOTAL_STATE; n++){
			pre_alpha_sum = 0.0;
			for (int n_pre = 0; n_pre < TOTAL_STATE; n_pre++){
				pre_alpha_sum += temp_paras->alpha[t-1][n_pre]*(hmm->transition[n_pre][n]);
			}
			temp_paras->alpha[t][n] = pre_alpha_sum * (hmm->observation[sequence[t]-65][n]);
		}
	}
    // beta t: T n: N
    for (int n = 0; n < TOTAL_STATE; n++){
		temp_paras->beta[sequence_len-1][n] = 1;
	}
	double post_beta_sum;
	for (int t = sequence_len-2; t>=0; t--){
		for (int n = 0; n < TOTAL_STATE; n++){
			post_beta_sum = 0.0;
			for (int n_post=0; n_post<TOTAL_STATE; n_post++){
				post_beta_sum += hmm->transition[n][n_post]*hmm->observation[sequence[t+1]-65][n_post]*temp_paras->beta[t+1][n_post];
			}
			temp_paras->beta[t][n] = post_beta_sum;
		}
	}

    // gamma t: T n: N
    double col_sum;
	for (int t=0; t<sequence_len; t++){
		col_sum=0.0;
		for (int n=0; n<TOTAL_STATE; n++){
			col_sum += temp_paras->alpha[t][n]*temp_paras->beta[t][n];
		}
		for (int n=0; n<TOTAL_STATE; n++){
			temp_paras->gamma[t][n] = temp_paras->alpha[t][n] * temp_paras->beta[t][n] / col_sum;
		}
	}

    // epsilon t: T n: N n: N
	double sum = 0.0;
    for (int t=0; t<sequence_len-1; t++){
		sum = 0.0;
		for (int n_i=0; n_i<TOTAL_STATE; n_i++){
			for (int n_j=0; n_j<TOTAL_STATE; n_j++){
				sum += (temp_paras->alpha[t][n_i]) * (hmm->transition[n_i][n_j]) * (hmm->observation[sequence[t+1]-65][n_j]) * (temp_paras->beta[t+1][n_j]);
			}
		}
		for (int n_i=0; n_i<TOTAL_STATE; n_i++){
			for (int n_j=0; n_j<TOTAL_STATE; n_j++){
				temp_paras->epsilon[t][n_i][n_j] = (temp_paras->alpha[t][n_i]) * (hmm->transition[n_i][n_j]) * (hmm->observation[sequence[t+1]-65][n_j]) * (temp_paras->beta[t+1][n_j]) / sum;
			}
		}
	}

    //--------------------------------sum!--------------------------------//
    sum_paras->times ++;
    
    //pi
    for (int n=0; n<TOTAL_STATE; n++){
		sum_paras->pi[n] += temp_paras->gamma[0][n];
    }

    //a_u
    double a_u_sum;
    for (int n_i=0; n_i<TOTAL_STATE; n_i++){
		for (int n_j=0; n_j<TOTAL_STATE; n_j++){
            a_u_sum = 0.0;
			for (int t=0; t<sequence_len-1; t++)
				a_u_sum += temp_paras->epsilon[t][n_i][n_j];
			sum_paras->a_u[n_i][n_j] += a_u_sum;
		}
	}

    //a_l
    double a_l_sum;
    for (int n=0; n<TOTAL_STATE; n++){
        a_l_sum = 0.0;
        for (int t=0; t<sequence_len-1; t++)
            a_l_sum += temp_paras->gamma[t][n];
        sum_paras->a_l[n] += a_l_sum;
	}
    
    //b_u
    double b_u_sum;
    for (int n=0; n<TOTAL_STATE; n++){
        for (int k=0; k<TOTAL_OBV; k++){
            b_u_sum = 0.0;
            for (int t=0; t<sequence_len; t++){
				if ((sequence[t]-65)==k){
					b_u_sum += temp_paras->gamma[t][n];
                }
            }
            sum_paras->b_u[n][k] += b_u_sum;
        }
    }

    //b_l
    double b_l_sum;
    for (int n=0; n<TOTAL_STATE; n++){
        b_l_sum = 0.0;
        for (int t=0; t<sequence_len; t++){
			b_l_sum += temp_paras->gamma[t][n];
        }
        sum_paras->b_l[n] += b_l_sum;
    }
    //--------------------------------------------------------------------//
}

void update_hmm(HMM *hmm, Sum_Paras *sum_paras){
    //pi
	for (int n=0; n<TOTAL_STATE; n++){
		hmm->initial[n] = sum_paras->pi[n]/sum_paras->times;
	}

    //a
	for (int n_i=0; n_i<TOTAL_STATE; n_i++){
		for (int n_j=0; n_j<TOTAL_STATE; n_j++){
			hmm->transition[n_i][n_j] = (sum_paras->a_u[n_i][n_j]) / (sum_paras->a_l[n_i]);
		}
	}

    //b
	for (int n_i=0; n_i<TOTAL_STATE; n_i++){
		for (int k=0; k<TOTAL_OBV; k++){
			hmm->observation[k][n_i] = (sum_paras->b_u[n_i][k])/(sum_paras->b_l[n_i]);
		}
	}
}
