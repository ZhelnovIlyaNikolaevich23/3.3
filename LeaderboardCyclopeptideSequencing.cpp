#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <ctime>

using std::map;
using std::cin;
using std::cout;
using std::map;
using std::string;
using std::vector;
using std::stoi;

int score_lin_pept(vector<int> &inp, vector<int> &spectrum);

struct score_peptides
{
	int score;
	vector<int> int_peptid;
	int mass;
	score_peptides()
	{
		mass = 0;
		score = 0;
	}

	score_peptides(const score_peptides& cop_pept)
	{
		mass = cop_pept.mass;
		score = cop_pept.score;
		int_peptid = cop_pept.int_peptid;
	}

	bool operator<(score_peptides &right)
	{
		return score < right.score;
	}
	bool operator<=(score_peptides &right)
	{
		return score <= right.score;
	}
	bool operator>(score_peptides &right)
	{
		return score > right.score;
	}
	bool operator>=(score_peptides &right)
	{
		return score >= right.score;
	}
	bool operator==(score_peptides &right)
	{
		return score == right.score;
	}
	bool operator!=(score_peptides &right)
	{
		return score != right.score;
	}
};


vector<int> pat = { 57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186 };

int score_lin_pept(vector<int> &inp, vector<int> &spectrum)
{
	vector<int> lin_pept = { 0 };
	int score = 0;
	for (int i = 0; i < inp.size(); i++)
	{
		for (int j = 0; j < inp.size() - i; j++)
		{
			int sum = 0;
			for (int k = j; k <= i + j; k++)
				sum = sum + inp[k];
			lin_pept.push_back(sum);
		}
	}
	map<int, int> kol_elem;
	for (int i = 0; i < spectrum.size(); i++)
		if (kol_elem.find(spectrum[i]) == kol_elem.end())
			kol_elem.insert({ spectrum[i], 1 });
		else
			kol_elem.at(spectrum[i])++;
	for (int i = 0; i < lin_pept.size(); i++)
	{
		if (kol_elem.find(lin_pept[i]) != kol_elem.end())
		{
			kol_elem.at(lin_pept[i])--;
			if (kol_elem.at(lin_pept[i]) >= 0)
				score++;
		}
	}
	return score;
}

vector <int> LeaderboardCyclopeptideSequencing(int N, vector<int> spectrum)
{
	score_peptides pept;
	vector<score_peptides> leaderboard;
	score_peptides out;
	leaderboard.push_back(pept);
	while (leaderboard.size())
	{
		vector<score_peptides> new_leaderboard;
		for (int i = 0; i < leaderboard.size(); i++)
		{
			for (int j = 0; j < pat.size(); j++)
			{
				score_peptides tmp;
				tmp.int_peptid = leaderboard[i].int_peptid;
				tmp.mass = leaderboard[i].mass;
				tmp.int_peptid.push_back(pat[j]);
				tmp.mass += pat[j];
				if (tmp.mass <= spectrum[spectrum.size() - 1])
				{
					tmp.score = score_lin_pept(tmp.int_peptid, spectrum);
					new_leaderboard.push_back(tmp);
				}
			}
		}
		std::sort(new_leaderboard.begin(), new_leaderboard.end());
		std::reverse(new_leaderboard.begin(), new_leaderboard.end());
		int i;
		for (i = N; (i < new_leaderboard.size() && new_leaderboard[i] == new_leaderboard[i - 1]); i++);
		if (new_leaderboard.size() <= N)
			leaderboard = new_leaderboard;
		else
		{
			leaderboard.resize(i);
			for (int k = 0; k < i; k++)
				leaderboard[k] = new_leaderboard[k];
			if (out.score < leaderboard[0].score)
			{
				out.score = leaderboard[0].score;
				out.mass = leaderboard[0].mass;
				out.int_peptid = leaderboard[0].int_peptid;
			}
		}
	}
	return out.int_peptid;
}



void main()
{
	unsigned int start_time = clock();
	int n;
	vector<int> plyx;
	cin >> n;
	char b;
	cin >> b;
	char c = '\0';
	int tmp = 0;
	while (c != '\n')
	{
		cin.get(c);
		if (c >= '0' && c <= '9')
		{
			tmp = tmp * 10 + (int)c - (int) '0';
		}
		else
		{
			plyx.push_back(tmp);
			tmp = 0;
		}
	}
	vector<int> obl = LeaderboardCyclopeptideSequencing(n, plyx);
	for (int i = 0; i < obl.size() - 1; i++)
		cout << obl[i] << '-';
	cout << obl[obl.size() - 1];
	unsigned int end_time = clock();
	cout << "\nTime: " << end_time - start_time;
}
