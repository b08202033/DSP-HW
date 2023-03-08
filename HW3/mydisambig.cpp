#include <stdio.h>
#include <string.h>
#include "Ngram.h"
#include <map>
#include <vector>

#ifndef MAX_CANDIDATE_NUM
    #define MAX_CANDIDATE_NUM 10000
#endif

#ifndef LEN
    #define LEN 500
#endif

using namespace std;

FILE *outFp;
VocabIndex UNK, a, b;
Vocab LM_vocab;
VocabIndex unigram_context[1] = {Vocab_None};
Ngram LM(LM_vocab, 2);

string startstr = "<s> ";
string endstr = "</s>\n";
const char* starting = startstr.c_str();
const char* ending = endstr.c_str();

map<string, vector<char*>> mapping;
map<char*, VocabIndex> index_table;
map<char*, VocabIndex>::iterator it;
map<char*, double> uni_prob_table;
map<char*, double>::iterator it2;
int l;

void buildDelta(double** delta, string& key){
	int n = mapping[key].size();
	delta[0] = new double [n];
	for(int i=0; i<n; i++){
		// initialized by unigram prob.
		it2 = uni_prob_table.find(mapping[key][i]);
		if(it2!=uni_prob_table.end()){delta[0][i] = it2->second;}
		else{
			VocabIndex wid;
			it = index_table.find(mapping[key][i]);
			if(it!=index_table.end()){wid = it->second;}
			else{
				wid = LM_vocab.getIndex(mapping[key][i]);
				if(wid==Vocab_None){wid = UNK;}
				index_table[mapping[key][i]] = wid;
			}
			double prob = LM.wordProb(wid, unigram_context);
			delta[0][i] = prob;
			uni_prob_table[mapping[key][i]] = prob;
		}
	}
}

int Viterbi(double** delta, int** back, VocabString w[LEN]){
	string previ(w[0]);
	int pre_num_cand = mapping[previ].size();
	int curr_num_cand;

	for(int i=1; i<l; i++){
		curr_num_cand = mapping[w[i]].size();
		delta[i] = new double [curr_num_cand];
		back[i] = new int [curr_num_cand];

		double best_prob;
		unsigned int best_prev;
		int j = 0;
		while(j<curr_num_cand){
			best_prob = -INFINITY;
			best_prev = 0;
			char* current_word = mapping[w[i]][j];

			it = index_table.find(current_word);
			if(it!=index_table.end()){b = it->second;}
			else{
				b = LM_vocab.getIndex(current_word);
				if(b==Vocab_None){b = UNK;}
				index_table[current_word] = b;
			}
			
			for(int k=0; k<pre_num_cand; k++){
				char* pre_word = mapping[previ][k];
				// read the bigram prob
				it = index_table.find(pre_word);
				if(it != index_table.end()){a = it->second;}
				else{
					a = LM_vocab.getIndex(pre_word);
					if(a==Vocab_None){a = UNK;}
					index_table[pre_word] = a;
				}
				
				VocabIndex bigram_context[2] = {a, Vocab_None};
				double prob = LM.wordProb(b, bigram_context) + delta[i-1][k];
				if(prob > best_prob){best_prev = k; best_prob = prob;}
			}
			delta[i][j] = best_prob; back[i][j] = best_prev;
			j++;
		}
		previ = w[i]; pre_num_cand = curr_num_cand;
	}
	return curr_num_cand;
}

void write_output(int& k, double** delta, int** back, VocabString w[LEN]){
	int psi = 0;
	for(int i=0; i<k; i++){if(delta[l-1][i]>delta[l-1][psi]){psi = i;}}
	vector<char*> output;
	for(int i=l-1; i>=0; i--){
		output.push_back(mapping[w[i]][psi]);
		if(i){psi = back[i][psi];}
	}
	fprintf(outFp, "%s", starting);
	for(int i=l-1; i>=0; i--){fprintf(outFp, "%s ", output[i]);}
	fprintf(outFp, "%s", ending);
}

void decode(char* s){
    VocabString w[LEN];
    l = Vocab::parseWords(s, w, LEN);
	string key(w[0]);
	// new the needed table
	int** back;
	back = new int* [l];
	double** delta;
	delta = new double* [l];
	// initialize the delta table
	buildDelta(delta, key);
	// Viterbi
	int k = Viterbi(delta, back, w);
	write_output(k, delta, back, w);

	return;
}

int main(int argc , char **argv){
	// open the file
    File be_decoded_file(argv[1] , "r");
	File map_file(argv[2] , "r");
	File LM_file(argv[3] , "r");
	// build the mapping
	char* buff;
	while((buff=map_file.getline())){
		VocabString words[MAX_CANDIDATE_NUM];
		int a = Vocab::parseWords(buff, words, MAX_CANDIDATE_NUM);

		vector<char*> c;
		for(int i=1; i<a; i++){
            char* ptr = new char [3];
			strncpy(ptr, words[i], 3);
			c.push_back(ptr);
        }
		
		mapping[words[0]] = c;
	}
	map_file.close();

	LM.read(LM_file);
	LM_file.close();
	UNK = LM_vocab.getIndex(Vocab_Unknown);
	
	char* str;
	outFp = fopen(argv[4], "w");
	while(str=be_decoded_file.getline()){decode(str);}

	be_decoded_file.close();
	fclose(outFp);

	return 0;
}

