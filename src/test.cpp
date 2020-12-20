#include "hmm.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
#define TOTAL_STATE 6
#define MODEL_NUM 5


double viterbi(HMM *, string);

int main(int argc, char **argv){

	const char* models_list_path = argv[1];
	const char* seq_path = argv[2];
	const char* output_result_path = argv[3];

	HMM hmms[MODEL_NUM];
	load_models(models_list_path, hmms, MODEL_NUM);

	double get_prob = 0.0;
	double best_prob = 0.0;
	double best_model = 1;

	string sequence="";
	ifstream ifs(seq_path, ifstream::in);
	ofstream ofs(output_result_path, ofstream::out);
	while(getline(ifs, sequence)){
		get_prob = 0.0;
		best_prob = 0.0;
		best_model = 1;
		for (int model_idx=0;model_idx<MODEL_NUM;model_idx++){
			get_prob = viterbi(&hmms[model_idx], sequence);
			if (get_prob >best_prob){
				best_prob = get_prob;
				best_model = model_idx + 1;
			}
		}
		ofs << "model_0" << best_model << ".txt " << best_prob << endl;
	}
	ifs.close();	
	ofs.close();
	return 0;
}

double viterbi(HMM *hmm, string sequence){
	int sequence_len = sequence.size();
	double delta[sequence_len][TOTAL_STATE];

	//init
	for (int n=0; n<TOTAL_STATE; n++){
		delta[0][n] = hmm->initial[n] * hmm->observation[sequence[0]-65][n];
	}

	//recursion
	double now_col = 0.0;
	double max_col = 0.0;
	for (int t=1; t<sequence_len; t++){
		for (int n=0; n<TOTAL_STATE; n++){
			max_col = 0.0;
			for (int n_pre=0; n_pre<TOTAL_STATE; n_pre++){
				now_col = delta[t-1][n_pre] * hmm->transition[n_pre][n];
				if (now_col>max_col){
					max_col = now_col;
				}
			delta[t][n] = max_col * hmm->observation[sequence[t]-65][n];
			}
		}
	}

	//termination
	double output_prob = 0.0;
	double now_prob = 0.0;
	for (int n=0; n<TOTAL_STATE; n++){
		now_prob = delta[sequence_len-1][n];
		if (now_prob > output_prob){
			output_prob = now_prob;
		}
	}
	return output_prob;
}
