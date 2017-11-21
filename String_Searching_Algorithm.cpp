#include "StringSearching.h"
#include <algorithm>
#include <chrono>	
#include <iomanip>
#include <iostream>
#include <cassert>

const unsigned CHARS = 256;


bool StringSearching::equalAnswer(const vector<unsigned> & A1, const vector<unsigned> & A2) {
	return A1.size() == A2.size() && A1 == A2;
}

void StringSearching::naiveStringMatcher(const string& pat) {
	unsigned n = text.length();
	unsigned m = pat.length();

	for (unsigned s = 0; s <= n - m; s++) {
		unsigned k;
		for (k = 0; k < m; k++) {
			if (pat[k] != text[s + k])
				break;
		}
		if (k == m) {
			answerNaive.push_back(s);
		}
	}
}

unsigned exp_mod(unsigned x, unsigned n, unsigned m) {
	if (n == 0)
		return 1;
	else if (n == 1)
		return x % m;
	else {
		unsigned square = (x*x) % m;
		unsigned exp = exp_mod(square, n / 2, m);

		if (n % 2 == 0)
			return exp % m;
		else
			return (exp*x) % m;
	}
}

unsigned hash_2(const string & str, unsigned len, unsigned b, unsigned m) {
	unsigned value = 0;
	unsigned power = 1;
	for (int i = len - 1; i >= 0; i--) {
		value += (power * str[i]);
		value %= m;
		power *= b;
		power %= m;
	}
	return value;
}

void StringSearching::rabinKarpMatcher(const string& pat) {
	unsigned m = pat.length();
	unsigned n = text.length();
	//long long h = 1, p = 0, t = 0;
	unsigned d = alphabet.size();
	unsigned q = 257;
	unsigned b = 1024;
	unsigned h = exp_mod(q, m, b);
	unsigned p = hash_2(pat, pat.size(), q, b);
	unsigned t = hash_2(text, pat.size(), q, b);
	/*for (unsigned j = 1; j < m; j++) {
	h = (d*h) % q;
	}*/

	/*for (unsigned j = 0; j < m; j++) {
	p = (p*d + pat[j]) % q;
	t = (t*d + text[j]) % q;
	}*/

	for (unsigned i = 0; i <= n - m; i++) {
		if (p == t) {
			unsigned j = 0;
			for (; j < m; j++) {
				if (text[i + j] != pat[j])
					break;
			}
			if (j == m) {
				answerRKM.push_back(i);
				//cout << i << endl;
			}
		}

		if (i < n - m) {
			t *= q;
			t -= ((text[i] * h) % b);
			t += text[i + m];
			t %= b;
			//t = (d*(t - text[i] * h) + text[i + m]) % q;
			if (t < 0) {
				t += b;
			}
		}
	}

}

//void StringSearching::rabinKarpMatcher(const string& pat) {
//	unsigned m = pat.length();
//	unsigned n = text.length();
//	long long h = 1, p = 0, t = 0;
//	unsigned d = alphabet.size();
//	unsigned q = 3355439;
//
//	for (unsigned j = 1; j < m; j++) {
//		h = (d*h) % q;
//	}
//
//	for (unsigned j = 0; j < m; j++) {
//		p = (p*d + pat[j]) % q;
//		t = (t*d + text[j]) % q;
//	}
//
//	for (unsigned i = 0; i <= n - m; i++) {
//		if (p == t) {
//			unsigned j = 0;
//			for (; j < m; j++) {
//				if (text[i + j] != pat[j])
//					break;
//			}
//			if (j == m)
//				answerRKM.push_back(i);
//		}
//
//		if (i < n - m) {
//			t = (d*(t - text[i] * h) + text[i + m]) % q;
//			if (t < 0) {
//				t += q;
//			}
//		}
//	}
//}

unsigned  StringSearching::nextState(const string& pat, const unsigned m, const unsigned state, const unsigned x) {

	if (state < m && x == pat[state])
		return state + 1;

	for (unsigned ns = state; ns > 0; ns--) {
		if (pat[ns - 1] == x) {
			unsigned i = 0;
			for (; i < ns - 1; i++)
				if (pat[i] != pat[state - ns + 1 + i])
					break;
			if (i == ns - 1)
				return ns;
		}
	}

	return 0;
}

void  StringSearching::fillTable(const string& pat, const unsigned m, unsigned **table) {
	for (unsigned state = 0; state <= m; state++) {
		for (unsigned x = 0; x < alphabet.size(); x++) {
			unsigned j = static_cast<unsigned>(alphabet[x]);
			table[state][j] = nextState(pat, m, state, j);
		}
	}
}

void StringSearching::finiteAutomationMatcher(const string& pat) {
	unsigned m = pat.length();
	unsigned n = text.length();

	unsigned **table = new unsigned*[m + 1];
	for (unsigned i = 0; i < m + 1; i++)
		table[i] = new unsigned[CHARS];

	fillTable(pat, m, table);

	unsigned state = 0;
	for (unsigned i = 0; i < n; i++) {
		state = table[state][text[i]];
		if (state == m) {
			answerSA.push_back(i - m + 1);
		}
	}
	for (unsigned i = 0; i < m + 1; i++)
		delete[] table[i];
	delete[] table;
}

void StringSearching::computePrefixFunction(const string& pat, vector<unsigned> &pref_func) {
	unsigned m = pat.length();
	pref_func[0] = 0;
	unsigned k = 0;
	for (unsigned q = 1; q < m; q++) {
		while (k > 0 && pat[k] != pat[q])
			k = pref_func[k - 1];
		if (pat[k] == pat[q])
			k++;
		pref_func[q] = k;
	}
}

void StringSearching::kmpMatcher(const string& pat) {
	unsigned m = pat.length();
	unsigned n = text.length();
	vector<unsigned> pref_func(m);
	computePrefixFunction(pat, pref_func);
	unsigned q = 0;
	for (unsigned i = 0; i < n; i++) {
		while (q > 0 && pat[q] != text[i])
			q = pref_func[q-1];
		if (pat[q] == text[i])
			q++;
		if (q == m) {
			answerKMP.push_back(i - m + 1);
			q = pref_func[q-1];
		}
	}
}

//char randFill(vector<char> alphabet) {
//	return alphabet[rand() % alphabet.size()];
//}

void StringSearching::changeText(const unsigned length) {
	if (patterns.size() == 0) {
		num_patterns = 1 + rand() % max_num_patterns;
		len_pattern = 5;//+ rand() % length
	}
	else {
		for (auto && pattern : patterns) {
			pattern.clear();
		}
		patterns.clear();
	}
	patterns.resize(num_patterns);
	auto comp = [&]()-> char {
		char ch = alphabet[rand() % alphabet.size()];
		return ch;
	};
	for (auto && pattern : patterns) {
		pattern.resize(len_pattern);
		generate(pattern.begin(), pattern.end(), comp);
	}
	text.clear();
	text.resize(length);
	generate(text.begin(), text.end(), comp);
}

void StringSearching::changePattern(const unsigned length) {
	if (patterns.size() == 0) {
		len_text = length * 2 + rand() % max_len;
		num_patterns = 1 + rand() % max_num_patterns;
	}
	else {
		for (auto && pattern : patterns) {
			pattern.clear();
		}
		patterns.clear();
	}
	patterns.resize(num_patterns);
	auto comp = [&]()-> char {
		char ch = alphabet[rand() % alphabet.size()];
		return ch;
	};
	for (auto && pattern : patterns) {
		pattern.resize(length);
		generate(pattern.begin(), pattern.end(), comp);
	}

	text.clear();
	text.resize(len_text);
	generate(text.begin(), text.end(), comp);

}

void StringSearching::changeNumPatterns(const unsigned number) {
	if (patterns.size() == 0) {
		len_text = 20 + rand() % max_len;
		len_pattern = 5;// + rand() % len_text
	}
	else {
		for (auto && pattern : patterns) {
			pattern.clear();
		}
		patterns.clear();
	}

	patterns.resize(number);
	auto comp = [&]()-> char {
		char ch = alphabet[rand() % alphabet.size()];
		return ch;
	};
	for (auto && pattern : patterns) {
		pattern.resize(len_pattern);
		generate(pattern.begin(), pattern.end(), comp);
	}

	text.clear();
	text.resize(len_text);
	generate(text.begin(), text.end(), comp);

}

void StringSearching::investigate(string fix_value, const unsigned length_start, const unsigned length_end, const unsigned step, const unsigned num_iter, ostream &out) {
	using namespace std::chrono;
	
		if (length_end > max_len)
			throw "length_end > max_len";

		transform(fix_value.begin(), fix_value.end(), fix_value.begin(), tolower);

		string ch;
		void(StringSearching::*change)(unsigned);
		if (fix_value.compare("text") == 0) {
			change = &StringSearching::changeText;
			ch = "n";
		}
		else if (fix_value.compare("pattern") == 0) {
			change = &StringSearching::changePattern;
			ch = "m";
		}
		else if (fix_value.compare("num_pattern") == 0) {
			change = &StringSearching::changeNumPatterns;
			ch = "np";
		}
		else
			throw "fix_value doesn't exist";

		out << "#" << setw(5) << " " << left << setw(20) << ch << left << setw(20) << "Naive" << left << setw(20) << "RK"
			<< left << setw(20) << "SA" << left << setw(20) << "KMP" << left << setw(20) << "AC" << endl;

		cout << "#" << setw(5) << " " << left << setw(20) << ch << left << setw(20) << "Naive" << left << setw(20) << "RK"
			<< left << setw(20) << "SA" << left << setw(20) << "KMP" << "AC" << endl;

		for (unsigned i = length_start; i <= length_end; i += step) {
			out << left << setw(6) << " " << left << setw(20) << i;
			cout << left << setw(6) << " " << left << setw(20) << i;

			double aver_time_naive = 0.0;
			double aver_time_karp = 0.0;
			double aver_time_stAuto = 0.0;
			double aver_time_kmp = 0.0;
			double aver_time_corasick = 0.0;

			for (unsigned j = 0; j < num_iter; j++) {
				(this->*change)(i);

				for (auto && pattern : patterns) {
					auto t_start = high_resolution_clock::now();
					naiveStringMatcher(pattern);
					auto t_end = high_resolution_clock::now();
					auto time = duration_cast<duration<double, milli>>(t_end - t_start);
					aver_time_naive += time.count();
					answerNaive.clear();
				}

				//answerNaive.clear();

				for (auto && pattern : patterns) {
					auto t_start = high_resolution_clock::now();
					rabinKarpMatcher(pattern);
					auto t_end = high_resolution_clock::now();
					auto time = duration_cast<duration<double, milli>>(t_end - t_start);
					aver_time_karp += time.count();
					answerRKM.clear();
				}

				//assert(equalAnswer(answerNaive, answerRKM));
				//answerRKM.clear();

				for (auto && pattern : patterns) {
					auto t_start = high_resolution_clock::now();
					finiteAutomationMatcher(pattern);
					auto t_end = high_resolution_clock::now();
					auto time = duration_cast<duration<double, milli>>(t_end - t_start);
					aver_time_stAuto += time.count();
					answerSA.clear();
				}
				

				//assert(equalAnswer(answerNaive, answerSA));
				//answerSA.clear();

				for (auto && pattern : patterns) {
					auto t_start = high_resolution_clock::now();
					kmpMatcher(pattern);
					auto t_end = high_resolution_clock::now();
					auto time = duration_cast<duration<double, milli>>(t_end - t_start);
					aver_time_kmp += time.count();
					answerKMP.clear();
				}

				//assert(equalAnswer(answerNaive, answerKMP));
				//answerKMP.clear();
				//answerNaive.clear();

				//t_start = high_resolution_clock::now();
				////ahoCorasick();
				//t_end = high_resolution_clock::now();
				//time = duration_cast<duration<double, milli>>(t_end - t_start);
				//aver_time_corasick += time.count();
			}

			aver_time_naive /= num_iter;
			aver_time_karp /= num_iter;
			aver_time_stAuto /= num_iter;
			aver_time_kmp /= num_iter;
			//aver_time_corasick /= num_iter;

			out << left << setw(20) << aver_time_naive << left << setw(20) << aver_time_karp
				<< left << setw(20) << aver_time_stAuto << left << setw(20) << aver_time_kmp
				<< aver_time_corasick << endl;

			cout << left << setw(20) << aver_time_naive << left << setw(20) << aver_time_karp
				 << left << setw(20) << aver_time_stAuto << left << setw(20) << aver_time_kmp
				 << aver_time_corasick << endl;
		}

		if (ch == "n") 
			out << endl << left << setw(10) << "#" << "m = " << left << setw(10) << len_pattern << "num_pat = " << left << setw(10) << num_patterns << "num_iter = " << num_iter;
		else if(ch == "m") 
			out << endl << left << setw(10) << "#" << "n = " << left << setw(10) << len_text <<  "num_pat = " << left << setw(10) << num_patterns << "num_iter = " << num_iter;
		else 
			out << endl << left << setw(10) << "#" << "m = " << left << setw(10) << len_pattern << "n = " << left << setw(16) << len_text << "num_iter = " << num_iter;
		
		if (ch == "n") 
			cout << endl << left << setw(10) << "#" << "m = " << left << setw(10) << len_pattern << "num_pat = " << left << setw(10) << num_patterns << "num_iter = " << num_iter << endl << endl ;
		else if (ch == "m") 
			cout << endl << left << setw(10) << "#" << "n = " << left << setw(10) << len_text << "num_pat = " << left << setw(10) << num_patterns << "num_iter = " << num_iter << endl << endl;
		else 
		cout << endl << left << setw(10) << "#" << "m = " << left << setw(10) << len_pattern << "n = " << left << setw(16) << len_text << "num_iter = " << num_iter << endl << endl;

		for (auto && pattern : patterns) {
			pattern.clear();
		}
		patterns.clear();
		text.clear();
		len_pattern = 0;
		len_text = 0;

		answerNaive.clear();
		answerRKM.clear();
		answerSA.clear();
		answerKMP.clear();
		answerAC.clear();

}