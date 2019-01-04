#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <unistd.h>

using namespace std;

#define POP_SIZE 500
#define RANDOM_NUM ((float)rand() / (RAND_MAX + 1)) //returns a float between 0 & 1
#define CORSSOVER_RATE 0.9
#define MUTATION_RATE 0.1

float UPPER = 10.0f;
float LOWER = -10.0f;
float X[POP_SIZE * 4][4] = {};
float fit[POP_SIZE * 4] = {};
float fit_sum;
float globalFit[5] = {};
int decimal;
int extra_cross_cromosome = 0;
int extra_mutation_crocromosome = 0;

int getRange(); //empty
void getRandomX(int c[][35], int, int);
void cal_fitness();
void reproduction_rws(int c[][35], int c2[][35]);
void crossover(int c2[][35]);
void mutation(int c2[][35]);
void cal_x(int c2[][35]);
void sort(int c[][35], int c2[][35]);

void sort(int c[][35], int c2[][35]) //gets rid of extra cromosome
{
    for (int i = 0; i < POP_SIZE / 2; i++)
    {
        float tmp = fit[0];
        int tmpNum = 0;
        for (int j = 0; j < (POP_SIZE + extra_cross_cromosome + extra_mutation_crocromosome - 1); j++)
        {
            if (tmp > fit[j + 1])
            {
                tmp = fit[j + 1];
                tmpNum = j + 1;
            }
        }
        for (int m = 0; m < 32; m++)
        {
            c[i][m] = c2[tmpNum][m];
        }

        fit[tmpNum] = 0.0f;
    }
}

void cal_x(int c2[][35])
{
    for (int i = 0; i < (POP_SIZE + extra_cross_cromosome + extra_mutation_crocromosome); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            int multi_count = 1;
            int k = (8 * j + 8 - 1);
            float zxc = 0.0f;
            for (k; k >= 8 * j; k--)
            {
                zxc += c2[i][k] * multi_count;
                multi_count *= 10;
            }
            float zzz;
            zzz = (zxc - 10000000) / 1000000;
            X[i][j] = zzz;
            //       cout << "asd = " << X[i][j] << endl;
        }
        //    cout << "------" << endl;
    }
}

void mutation(int c2[][35])
{
    extra_mutation_crocromosome = 0; //extra gene initialization

    for (int i = 0; i < POP_SIZE + extra_cross_cromosome; i++)
    {
        int if_mutated = 0;

        for (int k = 0; k < 32; k++)
        {
            if (RANDOM_NUM < MUTATION_RATE) //���ܲv
            {
                int rand_num = 9 * rand() / (RAND_MAX); //0~9
                if_mutated++;

                if (if_mutated > 1)
                {
                    //do
                }
                else
                {
                    for (int j = 0; j < 32; j++)
                    {
                        c2[POP_SIZE + extra_cross_cromosome + extra_mutation_crocromosome][j] = c2[i][j];
                    }
                    c2[POP_SIZE + extra_cross_cromosome + extra_mutation_crocromosome][k] = rand_num;
                }
            }
        }
        if (if_mutated > 0)
        {
            extra_mutation_crocromosome++;
        }
    }
}

void crossover(int c2[][35])
{
    extra_cross_cromosome = 0;
    int mating_pool[POP_SIZE] = {};
    int count = 0;
    for (int i = 0; i < POP_SIZE; i++)
    {
        if (RANDOM_NUM <= CORSSOVER_RATE)
        {
            mating_pool[count] = i;
            count++;
        }
    }

    for (int j = 0; j < POP_SIZE; j += 2)
    {
        if (mating_pool[j + 1] == 0)
        {
            break;
        }
        else
        {
            float linear_prob = RANDOM_NUM;
            int mat1 = mating_pool[j], mat2 = mating_pool[j + 1];
            int rand_num = 31 * rand() / (RAND_MAX); //0 ~ 31
            for (int k = 0; k < rand_num; k++)
            {
                c2[POP_SIZE + extra_cross_cromosome][k] = c2[mat1][k] * linear_prob + (1 - linear_prob) * c2[mat2][k];
                c2[POP_SIZE + extra_cross_cromosome + 1][k] = c2[mat1][k] * (1 - linear_prob) + c2[mat2][k] * linear_prob;
            }
            for (int m = rand_num; m < 32; m++)
            {
                c2[POP_SIZE + extra_cross_cromosome][m] = c2[mat2][m] * linear_prob + (1 - linear_prob) * c2[mat2][m];
                c2[POP_SIZE + extra_cross_cromosome + 1][m] = c2[mat1][m] * (1 - linear_prob) + c2[mat2][m] * linear_prob;
            }
            extra_cross_cromosome += 2;
        }
    }
}

void reproduction_rws(int c[][35], int c2[][35])
{
    float prob[POP_SIZE] = {};
    float trans_prob[POP_SIZE] = {};
    float trans_sum = 0.0;

    for (int i = 1; i < POP_SIZE; i++)
    {
        prob[i] = fit_sum / fit[i];
        trans_sum += prob[i];
    }
    trans_prob[0] = prob[0] / trans_sum;

    for (int i = 1; i < POP_SIZE; i++)
    {
        trans_prob[i] = prob[i] / trans_sum + trans_prob[i - 1];
    }

    for (int i = 0; i < POP_SIZE; i++)
    {
        float random = RANDOM_NUM;
        for (int j = 0; j < POP_SIZE; j++)
        {
            if (random <= trans_prob[j])
            {
                for (int k = 0; k < 32; k++)
                {
                    c2[i][k] = c[j][k];
                }
                c2[i][34] = (j + 1);
                //            cout << "repro = " << c2[i][34] << endl;
                break;
            }
        }
    }
}

void cal_fitness() //fx = 100*(x2 - x1^2)^2 || (1-x1)^2 +90(x4-x3^2)^2 || (1-x3)^2 + 10.1( (x2 - 1)^2 + (x4 - 1)^2 ) || 19.8*(x2-1)*(x4-1)
{
    fit_sum = 0;
    for (int i = 0; i < POP_SIZE + extra_cross_cromosome + extra_mutation_crocromosome; i++)
    {
        fit[i] = X[i][1] - pow(X[i][0], 2);
        fit[i] = 100 * pow(fit[i], 2);
        fit[i] += pow(1 - X[i][0], 2) + 90 * pow(X[i][3] - X[i][2] * X[i][2], 2);
        fit[i] += pow(1 - X[i][2], 2) + 10.1 * (pow(X[i][1] - 1, 2) + pow(X[i][3] - 1, 2));
        fit[i] += 19.8f * (X[i][1] - 1) * (X[i][3] - 1);
        //       cout << "i = " << i + 1 << " -> " << fit[i] << endl;
        fit_sum += fit[i];
        if (fit[i] < globalFit[0])
        {
            globalFit[0] = fit[i];
            globalFit[1] = X[i][0];
            globalFit[2] = X[i][1];
            globalFit[3] = X[i][2];
            globalFit[4] = X[i][3];
        }
    }
}

void getRandomX(int c[][35], int addCroNum, int start) //�٥�����Xrange //addcroNum = �X��cross
{
    int count = start;
    while (count < addCroNum)
    {
        //      cout << "count = " << count + 1 << endl;
        for (int i = 0; i < 4; i++)
        {
            float x = (rand() * (UPPER - LOWER) / RAND_MAX + LOWER);
            //         cout << setprecision(8) << "x = " << x << endl;
            x *= 1000000;
            x += 10000000;
            int xToInt = (int)x;
            float clone = xToInt;
            X[count][i] = (clone - 10000000) / 1000000;
            //         cout << setprecision(8) << "X[count][i] = " << X[count][i] << endl;
            int j = (8 * i + 8 - 1);
            for (j; j >= 8 * i; j--)
            {
                c[count][j] = xToInt % 10;
                xToInt /= 10;
                //      cout << c[count][j];
            }
            //           cout << endl;
        }

        //       cout << "-----" << endl;
        count++;
    }
}

float get_rand(float lower, float upper)
{
    return rand() * (upper - lower) / RAND_MAX + lower;
}

int main()
{
    globalFit[0] = 10000;
    int c[POP_SIZE * 3][35] = {};
    int c2[POP_SIZE * 3][35] = {};
    int let;
    srand((unsigned)time(NULL) + getpid());

    cin >> let;
    int cao = let;
    getRandomX(c, POP_SIZE, 0);
    cal_fitness();
    //while
    while (let > 0)
    {
        reproduction_rws(c, c2);
        crossover(c2);
        mutation(c2);
        cal_x(c2);
        cal_fitness();
        sort(c, c2);
        getRandomX(c, POP_SIZE, POP_SIZE / 2);
        extra_cross_cromosome = 0;
        extra_mutation_crocromosome = 0;
        let--;
        cout << "Iteration " << cao - let << " : " << globalFit[0] << endl;
    }
    //  cout << "fit_sum = " << fit_sum << endl;
    cout << endl;
    cout << "global = " << globalFit[0] << endl;
    cout << setprecision(7) << "X = " << globalFit[1] << " , " << globalFit[2] << " , " << globalFit[3] << " , " << globalFit[4] << endl;

    system("PAUSE");
}
//http://it-life.puckwang.com/2016/03/ccrandsrand.html
//http://it-easy.tw/cout-float/
//cout << setprecision(8) << get_rand(LOWER, UPPER) << endl;