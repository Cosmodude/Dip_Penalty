#include<iostream>
#include<fstream>
#include<iomanip>
#include <string>
#include<cstring>
#include <random>
#include <ctime>
#include <vector>

#define PI 3.14159265

using namespace std;
#pragma comment(linker, "/STACK:10000000")

class Galaxy
{
public:
	double ra;
	double dec;
	long long cw, ccw;
	double cosine(double alpha, double dzeta)
	{
		double cosine = cos(dzeta * PI / 180) * cos(dec * PI / 180) * cos((alpha - ra) * PI / 180) + sin(dzeta * PI / 180) * sin(dec * PI / 180);
		return cosine;
	}
};


inline void cnt_lines(string file,int c)
{
	ifstream fp(file.c_str());
	int g = 0;
	string str;
	if (!fp.is_open())
	{
		cout << "not opened" << endl;
	}
	while (fp >> g)
	{
		c++;
		getline(fp, str);
	}
}

inline void rand_one(int c, vector <vector<int>>(&mass))
{
	mass.clear();
	vector <int> vec;
	int g = 0;

	//cout << c << endl;
	for (int i = 0; i < 1001; i++)
	{
		vec.clear();
		for (int j = 0; j < c; j++)
		{
			srand(static_cast<unsigned int>(time(0)));
			int val = rand();
			if (val % 2 == 1)
			{
				g = 1;
			}
			else
			{
				g = -1;
			}
			vec.push_back(g);
		}
		mass.push_back(vec);
	}
}

inline void read_t(vector <Galaxy>(&mass), string file)
{
	srand(static_cast<unsigned int>(time(0)));

	mass.clear();
	string str;
	Galaxy g;
	ifstream fp(file.c_str());
	if (!fp.is_open())
	{
		cout << "not opened" << endl;
	}
	while (fp >> g.cw)
	{
		fp.get();
		fp >> g.ccw;
		fp.get();
		fp >> g.ra;
		fp.get();
		fp >> g.dec;
		getline(fp, str);
		mass.push_back(g);
	}
}

double func_real(double alpha, double dzeta, vector <Galaxy>(&mass), vector<vector<int>>(&one))
{
	double sum = 0;
	vector<double> cosine;

	for (int i = 0; i < mass.size(); i++)
	{
		cosine.push_back(mass[i].cosine(alpha, dzeta));
//random data use
		if (abs(cosine[i]) * (one[1000][i]) == cosine[i])
		{
			sum++;
		}
	}
	//cout << sum << endl;
	return pow(mass.size() - sum, 2) / mass.size();
}

double func_rand(double alpha, double dzeta, vector<Galaxy>(&mass), vector <vector<int>>(&one), double(&sumr))
{
	double cra[1000];
	double sigma = 0;
	double sum_cra = 0;
	vector<double> cosine;
	int count = 0;
	for (int i = 0; i < 1000; i++)
	{
		cosine.clear();
		count = 0;
		for (int j = 0; j < mass.size(); j++)
		{
			cosine.push_back(mass[j].cosine(alpha, dzeta));
			//random data use
			if (abs(cosine[j]) * (one[i][j]) == cosine[j])
			{
				count++;
			}
		}
		cra[i] = pow(mass.size() - count, 2) / mass.size();
		sum_cra += cra[i];
	}
	//cout << cra[2]<<endl;
	sumr = sum_cra / 1000;

	double x = 0;
	for (int i = 0; i < 1000; i++)
	{
		x += pow(sumr - cra[i], 2);
	}

	sigma = sqrt(x / 1000);

	return sigma;
}

int main()
{
	srand(static_cast<unsigned int>(time(0)));

	double n=0,m=0,cw=0,tr=0;
	vector <Galaxy> tdata;
	read_t(tdata, "catalog.csv");
	vector <vector <int>> fdata;
	ofstream res("results.txt");
	ofstream pic("picture.txt");
	
	while (tr < 101) 
	{
		double sumr = 0,pr=0, temp=0;
		res << "in process " << n+1 << endl;
		rand_one( tdata.size(),fdata );
		ofstream fout("pen_dip.txt");
		for (double ra = 0; ra < 361; ra += 5)
		{
			for (double dec = -90; dec < 91; dec += 5)
			{
				temp=0;
				temp = func_rand(ra, dec, tdata, fdata, sumr);
				pr = (abs(sumr - func_real(ra, dec, tdata,fdata)) / temp);
				fout << ra << "\t" << dec << "\t" << pr << endl;
				if(n==0){ pic << ra << "\t" << dec << "\t" << pr << endl; }
				m++;
				if (pr > 2.4)
				{
					tr = 200;
					res << "ra = " << ra << "  " << "dec = " << dec << "  "<< "pr = "<<pr<< endl;
				}
			}
			fout << endl;
			if (n == 0) { pic << endl; }
		}
		
		for (int i = 0; i < tdata.size()-1; i++)
			{
			if (fdata[1000][i] == 1) {cw++;}
			}
	
		res <<"cw = "<< cw << endl;
		cw = 0;
		n++;
	}

	res<< endl <<"n "<< n <<endl;
	res << "m "<< m << endl;
	res << "done" << endl;

	system("pause");
	return 0;
}