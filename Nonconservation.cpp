#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

//����
double gamma = 1.4;
double CFL = 0.5;

//�ռ�
const int N = 31;
double dx = 0.1;
double x[N];
double A[N];

//ʱ��
int timestep = 1400;
double dt[N];
double dt_min;

//������
double rho[N], T[N], V[N];//��ʼֵ
double p[N],Ma[N],U2[N];//ѹ�������������������
double rho_t[N], T_t[N], V_t[N];//Ԥ�ⲽ
double rho_bar[N], T_bar[N], V_bar[N];//Ԥ��ֵ
double rho_t_bar[N], T_t_bar[N], V_t_bar[N];//������
double rho_t_av[N], T_t_av[N], V_t_av[N];//ƽ��ֵ
double rho2[N], T2[N], V2[N];//����ֵ


int main()
{
	int i, j;

	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;

	//��ʼ��
	for (i = 0; i < N ; i++)
	{
		x[i] = i * dx;
		A[i] = 1 + 2.2 * (x[i] - 1.5) * (x[i] - 1.5);
		rho[i] = 1 - 0.314 * x[i];
		T[i] = 1 - 0.2314 * x[i];
		V[i] = (0.1 + 1.09 * x[i]) * sqrt(T[i]);
	}


	file2.open("x=1.5m����������ʱ�䲽�ı仯.txt");
	file2 << "t\t" << "rho2\t" << "V2\t" << "T2\t" << "p\t" << "Ma\t" << "\n";
	file3.open("x=1.5m���в�����.txt");
	file3 << "t\t" << "rho_t_av\t" << "V_t_av\t"  << "\n";
	file4.open("dt_min.txt");
	for (j = 0; j < timestep; j++)
	{
		//ʱ�䲽���ļ���
		dt_min= CFL * dx / (sqrt(T[0]) + V[0]);
		for (i = 0; i < N; i++)
		{
			dt[i] = CFL * dx / (sqrt(T[i]) + V[i]);
			if (dt_min > dt[i])
			{
				dt_min = dt[i];
			}
		}

		//Ԥ�����ļ���
		for (i = 0; i < N - 1; i++)
		{
			rho_t[i] = (-V[i] * (rho[i + 1] - rho[i]) - rho[i] * (V[i + 1] - V[i]) - rho[i] * V[i] * (log(A[i + 1]) - log(A[i]))) / dx;
			V_t[i] = (-V[i] * (V[i + 1] - V[i]) - ((T[i + 1] - T[i]) + T[i] / rho[i] * (rho[i + 1] - rho[i])) / gamma) / dx;
			T_t[i] = (-V[i] * (T[i + 1] - T[i]) - (gamma - 1) * T[i] * ((V[i + 1] - V[i]) + V[i] * (log(A[i + 1]) - log(A[i])))) / dx;
		}

		//Ԥ��ֵ�ļ���
		for (i = 0; i < N - 1; i++)
		{
			rho_bar[i] = rho[i] + rho_t[i] * dt_min;
			V_bar[i] = V[i] + V_t[i] * dt_min;
			T_bar[i] = T[i] + T_t[i] * dt_min;
		}

		//�������ļ���
		for (i = 1; i < N - 1; i++)
		{
			rho_t_bar[i] = (-V_bar[i] * (rho_bar[i] - rho_bar[i - 1]) - rho_bar[i] * (V_bar[i] - V_bar[i - 1]) - rho_bar[i] * V_bar[i] * (log(A[i]) - log(A[i - 1]))) / dx;
			V_t_bar[i] = (-V_bar[i] * (V_bar[i] - V_bar[i - 1]) - ((T_bar[i] - T_bar[i - 1]) + T_bar[i] / rho_bar[i] * (rho_bar[i] - rho_bar[i - 1])) / gamma) / dx;
			T_t_bar[i] = (-V_bar[i] * (T_bar[i] - T_bar[i - 1]) - (gamma - 1) * T_bar[i] * ((V_bar[i] - V_bar[i - 1]) + V_bar[i] * (log(A[i]) - log(A[i - 1])))) / dx;
		}

		//ƽ��ֵ�ļ���
		for (i = 1; i < N - 1; i++)
		{
			rho_t_av[i] = (rho_t[i] + rho_t_bar[i]) / 2;
			V_t_av[i] = (V_t[i] + V_t_bar[i]) / 2;
			T_t_av[i] = (T_t[i] + T_t_bar[i]) / 2;
		}

		//����ֵ�ļ���
		for (i = 1; i < N - 1; i++)
		{
			rho2[i] = rho[i] + rho_t_av[i] * dt_min;
			V2[i] = V[i] + V_t_av[i] * dt_min;
			T2[i] = T[i] + T_t_av[i] * dt_min;
		}


		//�߽���ֵ
		V2[0] = 2 * V2[1] - V2[2];
		rho2[0] = 1;
		T2[0] = 1;

		V2[N - 1] = 2 * V2[N - 2] - V2[N - 3];
		rho2[N - 1] = 2 * rho2[N - 2] - rho2[N - 3];
		T2[N - 1] = 2 * T2[N - 2] - T2[N - 3];

		for (i = 0; i < N; i++)
		{
			V[i] = V2[i];
			rho[i] = rho2[i];
			T[i] = T2[i];
		}

		//ѹ���������
		for (i = 0; i < N ; i++)
		{
			p[i] = rho2[i] * T2[i];
			Ma[i] = V2[i] / sqrt(T2[i]);
			U2[i] = rho2[i] * A[i] * V2[i];
		}
		
		file2 << j <<"\t"<< rho2[15] << "\t" << V2[15] << "\t" << T2[15] << "\t" << p[15] << "\t" << Ma[15] << endl;
		file3 << j << "\t" << abs(rho_t_av[15]) << "\t" << abs(V_t_av[15]) << "\t" << endl;
		file4 << j << "\t" << dt_min << endl;
	}
	file2.close();
	file3.close();
	file4.close();

	file1.open("1400ʱ�䲽�����������.txt");
	file1 << "x\t" << "rho2\t" << "V2\t" << "T2\t" <<"p\t" << "Ma\t" << "U2\t" "\n";
	for (i = 0; i < N; i++)
	{
		file1 << x[i] << "\t" << rho2[i] << "\t" << V2[i] << "\t" << T2[i] << "\t"<<p[i]<<"\t" << Ma[i] << "\t" << U2[i] <<endl;
	}
	file1.close();

}
