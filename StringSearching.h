#ifndef STRING_SEARCHING_H
#define STRING_SEARCHING_H

#include <string>
#include <vector>
#include <fstream>
using namespace std;

class StringSearching {
private:
	vector<char> alphabet;
	string text;
	vector<string> patterns;
	unsigned num_patterns;
	unsigned max_len;
	unsigned max_num_patterns;
	unsigned len_pattern;
	unsigned len_text;

	vector<unsigned> answerNaive;
	vector<unsigned> answerRKM;
	vector<unsigned> answerSA;
	vector<unsigned> answerKMP;
	vector<unsigned> answerAC;

	void changeText(unsigned length);
	void changePattern(unsigned length);
	void changeNumPatterns(unsigned number);
	void fillTable(const string& pat, unsigned m, unsigned **table);
	unsigned nextState(const string& pat, unsigned m, unsigned state, unsigned x);
	void computePrefixFunction(const string& pat, vector<unsigned> &pref_func);
	bool equalAnswer(const vector<unsigned> & A1, const vector<unsigned> & A2);
public:
	StringSearching(vector<char> alphabet, unsigned max_len, unsigned max_num_patterns) {
		this->alphabet = alphabet;
		this->max_len = max_len;
		this->max_num_patterns = max_num_patterns;
	}

	void investigate(string fix_value, const unsigned length_start, const unsigned length_end, const unsigned step, const unsigned num_iter, ostream &out);

	void naiveStringMatcher(const string& pat);
	void rabinKarpMatcher(const string& pat);
	void finiteAutomationMatcher(const string& pat);
	void kmpMatcher(const string& pat);


};
#endif // STRING_SEARCHING_H
