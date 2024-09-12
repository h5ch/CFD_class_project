#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

//常量
double gamma = 1.4;
double CFL = 0.5;

//空间
const int N = 31;
double dx = 0.1;
double x[N];
double A[N];

//时间
int timestep = 1400;
double dt[N];
double dt_min;

//物理量
double rho[N], T[N], V[N];//初始值
double p[N],Ma[N],U2[N];//压力，马赫数，质量流量
double rho_t[N], T_t[N], V_t[N];//预测步
double rho_bar[N], T_bar[N], V_bar[N];//预测值
double rho_t_bar[N], T_t_bar[N], V_t_bar[N];//修正步
double rho_t_av[N], T_t_av[N], V_t_av[N];//平均值
double rho2[N], T2[N], V2[N];//修正值


int main()
{
	int i, j;

	ofstream file1;
	ofstream file2;
	ofstream file3;
	ofstream file4;

	//初始化
	for (i = 0; i < N ; i++)
	{
		x[i] = i * dx;
		A[i] = 1 + 2.2 * (x[i] - 1.5) * (x[i] - 1.5);
		rho[i] = 1 - 0.314 * x[i];
		T[i] = 1 - 0.2314 * x[i];
		V[i] = (0.1 + 1.09 * x[i]) * sqrt(T[i]);
	}


	file2.open("x=1.5m处物理量随时间步的变化.txt");
	file2 << "t\t" << "rho2\t" << "V2\t" << "T2\t" << "p\t" << "Ma\t" << "\n";
	file3.open("x=1.5m处残差曲线.txt");
	file3 << "t\t" << "rho_t_av\t" << "V_t_av\t"  << "\n";
	file4.open("dt_min.txt");
	for (j = 0; j < timestep; j++)
	{
		//时间步长的计算
		dt_min= CFL * dx / (sqrt(T[0]) + V[0]);
		for (i = 0; i < N; i++)
		{
			dt[i] = CFL * dx / (sqrt(T[i]) + V[i]);
			if (dt_min > dt[i])
			{
				dt_min = dt[i];
			}
		}

		//预估步的计算
		for (i = 0; i < N - 1; i++)
		{
			rho_t[i] = (-V[i] * (rho[i + 1] - rho[i]) - rho[i] * (V[i + 1] - V[i]) - rho[i] * V[i] * (log(A[i + 1]) - log(A[i]))) / dx;
			V_t[i] = (-V[i] * (V[i + 1] - V[i]) - ((T[i + 1] - T[i]) + T[i] / rho[i] * (rho[i + 1] - rho[i])) / gamma) / dx;
			T_t[i] = (-V[i] * (T[i + 1] - T[i]) - (gamma - 1) * T[i] * ((V[i + 1] - V[i]) + V[i] * (log(A[i + 1]) - log(A[i])))) / dx;
		}

		//预估值的计算
		for (i = 0; i < N - 1; i++)
		{
			rho_bar[i] = rho[i] + rho_t[i] * dt_min;
			V_bar[i] = V[i] + V_t[i] * dt_min;
			T_bar[i] = T[i] + T_t[i] * dt_min;
		}

		//修正步的计算
		for (i = 1; i < N - 1; i++)
		{
			rho_t_bar[i] = (-V_bar[i] * (rho_bar[i] - rho_bar[i - 1]) - rho_bar[i] * (V_bar[i] - V_bar[i - 1]) - rho_bar[i] * V_bar[i] * (log(A[i]) - log(A[i - 1]))) / dx;
			V_t_bar[i] = (-V_bar[i] * (V_bar[i] - V_bar[i - 1]) - ((T_bar[i] - T_bar[i - 1]) + T_bar[i] / rho_bar[i] * (rho_bar[i] - rho_bar[i - 1])) / gamma) / dx;
			T_t_bar[i] = (-V_bar[i] * (T_bar[i] - T_bar[i - 1]) - (gamma - 1) * T_bar[i] * ((V_bar[i] - V_bar[i - 1]) + V_bar[i] * (log(A[i]) - log(A[i - 1])))) / dx;
		}

		//平均值的计算
		for (i = 1; i < N - 1; i++)
		{
			rho_t_av[i] = (rho_t[i] + rho_t_bar[i]) / 2;
			V_t_av[i] = (V_t[i] + V_t_bar[i]) / 2;
			T_t_av[i] = (T_t[i] + T_t_bar[i]) / 2;
		}

		//修正值的计算
		for (i = 1; i < N - 1; i++)
		{
			rho2[i] = rho[i] + rho_t_av[i] * dt_min;
			V2[i] = V[i] + V_t_av[i] * dt_min;
			T2[i] = T[i] + T_t_av[i] * dt_min;
		}


		//边界点的值
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

		//压力，马赫数
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

	file1.open("1400时间步后的流场参数.txt");
	file1 << "x\t" << "rho2\t" << "V2\t" << "T2\t" <<"p\t" << "Ma\t" << "U2\t" "\n";
	for (i = 0; i < N; i++)
	{
		file1 << x[i] << "\t" << rho2[i] << "\t" << V2[i] << "\t" << T2[i] << "\t"<<p[i]<<"\t" << Ma[i] << "\t" << U2[i] <<endl;
	}
	file1.close();

}
