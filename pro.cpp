#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
using namespace std;

void preKMP(string pattern, int f[])
{
    int m = pattern.length(), k;
    f[0] = -1;
    for (int i = 1; i < m; i++)
    {
        k = f[i - 1];
        while (k >= 0)
        {
            if (pattern[k] == pattern[i - 1])
                break;
            else
                k = f[k];
        }
        f[i] = k + 1;
    }
}

//check whether target string contains pattern

int KMP(string pattern, string target)
{
    int m = pattern.length();
    int n = target.length();
    int f[m];
    preKMP(pattern, f);
    int i = 0;
    int k = 0;
    while (i < n)
    {
        if (k == -1)
        {
            i++;
            k = 0;
        }
        else if (target[i] == pattern[k])
        {
            i++;
            k++;
            if (k == m)
                return i-m;
        }
        else
            k = f[k];
    }
    return -1;
}





const size_t alphabets = 52;

int align(const string &a, const string &b, int alpha_gap,
        int alpha[alphabets][alphabets], string &a_aligned,
        string &b_aligned);

void print2DVector(const vector<vector<int> > &A);

int min(int a, int b, int c);

int glbal(string a1,string b1)
{


    // Penalty for any alphabet matched with a gap
    int gap_penalty = 2;

    /*
     * alpha[i][j] = penalty for matching the ith alphabet with the
     *               jth alphabet.
     * Here: Penalty for matching an alphabet with anoter one is 1
     *       Penalty for matching an alphabet with itself is 0
     */
    int alpha[alphabets][alphabets];
    for (size_t i = 0; i < alphabets; ++i)
    {
        for (size_t j = 0; j < alphabets; ++j)
        {
            if (i == j&& i<26) alpha[i][j] = -3;
            else if(i == j&& i>=26)  alpha[i][j] = -50;
            else alpha[i][j] = 1;
        }
    }

    // Aligned sequences
    string a2, b2;
    int penalty = align(a1, b1, gap_penalty, alpha, a2, b2);
    cout << "Needleman-Wunsch Score: " << penalty << endl;
    //cout << "Aligned sequences: " << endl;
    //cout << a2 << endl;
    //cout << b2 << endl;
    return penalty;
}


int align(const string &a, const string &b, int alpha_gap,
        int alpha[alphabets][alphabets], string &a_aligned,
        string &b_aligned)
{
    size_t n = a.size();
    size_t m = b.size();

    vector<vector<int> > A(n + 1, vector<int>(m + 1));

    for (size_t i = 0; i <= m; ++i)
        A[0][i] = alpha_gap * i;
    for (size_t i = 0; i <= n; ++i)
        A[i][0] = alpha_gap * i;

    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];
            if((x_i >= 'a' && x_i <= 'z')&&(y_j >= 'a' && y_j <= 'z'))
            A[i][j] = min(A[i-1][j-1] + alpha[x_i - 'a'+ 26][y_j - 'a' + 26],
                          A[i-1][j] + alpha_gap,
                          A[i][j-1] + alpha_gap);
            else if((x_i >= 'a' && x_i <= 'z')||(y_j >= 'a' && y_j <= 'z'))
            A[i][j] = min(A[i-1][j-1] + 1,
                          A[i-1][j] + alpha_gap,
                          A[i][j-1] + alpha_gap);
            else
            A[i][j] = min(A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'],
                          A[i-1][j] + alpha_gap,
                          A[i][j-1] + alpha_gap);
        }
    }

    // print2DVector(A);

    a_aligned = "";
    b_aligned = "";
    size_t j = m;
    size_t i = n;
    for (; i >= 1 && j >= 1; --i)
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (x_i == y_j && A[i][j] == A[i-1][j-1] - 3)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
        if (x_i == y_j && A[i][j] == A[i-1][j-1] - 50)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
        else if (x_i != y_j && A[i][j] == A[i-1][j-1] + 1)
        {
            /*
             */
            a_aligned = x_i + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
        else if (A[i][j] == A[i-1][j] + alpha_gap)
        {
            a_aligned = x_i + a_aligned;
            b_aligned = '-' + b_aligned;
        }
        else
        {
            a_aligned = '-' + a_aligned;
            b_aligned = y_j + b_aligned;
            --j;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned = a[i-1] + a_aligned;
        b_aligned = '-' + b_aligned;
        --i;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned = '-' + a_aligned;
        b_aligned = b[j-1] + b_aligned;
        --j;
    }

    return A[n][m];
}

int min(int a, int b, int c)
{
    if (a <= b && a <= c)
        return a;
    else if (b <= a && b <= c)
        return b;
    else
        return c;
}


void print2DVector(const vector<vector<int> > &A)
{
    for (auto i : A)
    {
        for (auto j : i)
            cout << j << " ";
        cout << endl;
    }
}



int main()
{
    string str1,str2,min_seq;
    int no_of_core,s,p,no_of_seq,min_score = 0,penalty,min =0;
    cout << "Enter a sequence:" << endl;
    cin >> str1;
    fstream fp;
    fp.open("foldlib.txt",ios::in);
    fp >> no_of_seq;
  //  cout << no_of_seq<<"\n";
    for(int k = 1;k<=no_of_seq;k++){
        fp>>str2;
        string temp1 = str1;
        string temp2 = str2;
  //      cout<<str1<<"\n"<<str2<<"\n";
        fp>>no_of_core;
        string str;
        for(int i=1;i<=no_of_core;i++)
        {
            fp>>str;
            s=KMP(str,str1);
            if(s!=-1){
            for(int j=0;j<str.size();j++)
                temp1[j+s]=temp1[j+s] + 32;
            p=KMP(str,str2);
            for(int j=0;j<str.size();j++)
                temp2[j+p]=temp2[j+p] + 32;

            }
        }
        cout<< temp1 << "\n" << temp2 << "\n";
        penalty = glbal(temp1,temp2);
        if(min_score > penalty){
            min_score = penalty;
            min_seq = str2;
        }
    }
    cout << "\nMin Score:"<< min_score << endl << "Best Fit sequence:" << min_seq << endl;

return 0;
}
