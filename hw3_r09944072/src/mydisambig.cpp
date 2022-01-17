// man page: http://www.speech.sri.com/projects/srilm/manpages/
/*
todo:
1. build_map   ---> 他是一個VocabMap
2. readin language model ---> 他是一個N-gram
3. read text to Viterbi
*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "Ngram.h"
#define PROB_MIN -20
#define VERY_NEG -1000000000
#define MAX_LINE_LENGTH 200000

bool g_map_built = false;
volatile int g_ngram_order = 2;

class MyMapping {
    public:
        map<int, vector<string>> internal_map;

        MyMapping() {}
        ~MyMapping() {}
        int hash_string(string s) {
            if (s == "<s>") {
                return 0x11223344;
            }
            else if (s == "</s>") {
                return 0x55667788;
            }
            assert (s.length() == 2);
            return (int(s[0]) * 256) + int(s[1]);
        }

        vector<string> getValue(string s) {
            int key = hash_string(s);
            return internal_map[key];
        }

        void setValue(string s, vector<string> &vec) {
            int key = hash_string(s);
            internal_map[key] = vec;
        }

        void build_map(char *map_filename) {
            if (g_map_built) return;
            g_map_built = true;

            // Build ZhuYin-Big5 map
            FILE *mapfile = fopen(map_filename, "r");
            char raw_data[MAX_LINE_LENGTH] = {0};

            while (fgets(raw_data, MAX_LINE_LENGTH, mapfile)) {
                string data(raw_data);
                string idx = data.substr(0, 2);

                vector<string> words;
                words.push_back(data.substr(3, 2));
                int pos = 2;
                while (1) {
                    pos = data.find(" ", pos+1);
                    if (pos == string::npos) break;
                    words.push_back(data.substr(pos+1, 2));
                }
                setValue(idx, words);
            }

            // add begin and terminal tag
            vector<string> b (1, "<s>"), t (1, "</s>");
            setValue("<s>", b);
            setValue("</s>", t);
        }
};

Vocab voc;
Ngram lm(voc, 2);
MyMapping my_map;

double ngram_prob(string word, string prev, int ngram_order) {
    VocabIndex wid = voc.getIndex(word.c_str());
    if (wid == Vocab_None)
        return PROB_MIN;
    if (ngram_order == 2) {
        VocabIndex context[] = {voc.getIndex(prev.c_str()), Vocab_None};
        return lm.wordProb(wid, context);
    }
    else if (ngram_order == 3) {
        // TODO: implement trigram
        VocabIndex context[] = {voc.getIndex(prev.c_str()), voc.getIndex(prev.c_str()), Vocab_None};
        return lm.wordProb(wid, context);
    }
    return -1.0;
}

vector<string> Viterbi(char *argv[], vector<string> &words) {
    my_map.build_map(argv[2]);

    VocabIndex wid;
    VocabIndex context[] = {Vocab_None, Vocab_None};
    vector<vector<double>> v;
    vector<vector<int>> bt;

    // initialize
    for (int i = 0; i < words.size(); i++) {
        int size;
        if (i == 0) {
            size = 1;
        }
        else {
            vector<string> vec = my_map.getValue(words[i]);
            size = vec.size();
        }
        vector<double> w (size, PROB_MIN);
        vector<int> b (size, 0);
        v.push_back(w);
        bt.push_back(b);
    }
    v[0][0] = 1;

    // recursion
    for (int i = 1; i < words.size(); i++) {
        vector<string> vec = my_map.getValue(words[i]);
        vector<string> vec_prev = my_map.getValue(words[i-1]);

        int j = 0;
        for (auto it = vec.begin(); it != vec.end(); j++, it++) {
            double prob = VERY_NEG, max_prob = VERY_NEG;
            for (int k = v[i-1].size()-1; k >= 0; k--) {
                prob = ngram_prob(*it, vec_prev[k], g_ngram_order);
                if (v[i-1][k] + prob > max_prob) {
                    max_prob = v[i-1][k] + prob;
                    v [i][j] = max_prob;
                    bt[i][j] = k;
                }
            }
        }
    }

    // termination
    int idx = words.size()-1, which = -1;
    double max_prob = VERY_NEG;
    vector<string> vec = my_map.getValue(words[idx]);
    int w_size = vec.size();
    for (int i = w_size-1; i >= 0; i--) {
        if (v[idx][i] > max_prob) {
            max_prob = v[idx][i];
            which = i;
        }
    }
    if (which == -1) {
        printf("error: didn't find the best\n");
        exit(1);
    }

    // backtrace
    vector<string> ans (words.size(), "");
    int best[MAX_LINE_LENGTH] = {0};
    for (int i = idx; i >= 0; i--) {
        vector<string> vec = my_map.getValue(words[i]);
        best[i] = which;
        which = bt[i][which];
    }

    // set the best word
    for (int i = 0; i <= idx; i++) {
        vector<string> vec = my_map.getValue(words[i]);
        ans[i] = vec[best[i]];
    }
    return ans;
}


int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: ./mydisambig text FILE-MAP LM output.txt\n");
        return 0;
    }

    File lmFile(argv[3], "r");
    lm.read(lmFile);
    lmFile.close();

    FILE *infile = fopen(argv[1], "r");
    char data[MAX_LINE_LENGTH] = {0};
    if (!infile) {
        printf("error reading file %s", argv[1]);
        return 0;
    }

    FILE *outfile = fopen(argv[4], "w");
    if (!outfile) {
        printf("error reading file %s", argv[4]);
        return 0;
    }

    while (fgets(data, MAX_LINE_LENGTH, infile)) {
        // extract words
        vector<string> words;
        words.push_back("<s>");
        for (int i = 0; i < MAX_LINE_LENGTH; i++) {
            if (data[i] == ' ') continue;
            if (data[i] == '\n' || data[i] == '\0') break;
            // 一個中文字：2byte
            char substr[3] = {0};
            substr[0] = data[i];
            substr[1] = data[i+1];
            string s(substr);
            words.push_back(s); // string substr (size_t pos = 0, size_t len = npos) const;
            i++;
        }
        words.push_back("</s>");

        // Viterbi algorithm
        vector<string> result = Viterbi(argv, words);

        result.pop_back();
        for (auto it = result.begin(); it != result.end(); it++) {
            fprintf(outfile, "%s ", (*it).c_str());
        }
        fputs("</s>\n", outfile);

    }

    // closing files
    if (fclose(infile)) {
        printf("error closing file %s", argv[1]);
    }
    if (fclose(outfile)) {
        printf("error closing file %s", argv[4]);
    }
    return 0;
}
