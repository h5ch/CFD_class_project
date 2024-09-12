#include<iostream>;
#include<cmath>;
#include<fstream>;
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
double rho[N], T[N], V[N];//原始变量
double U1[N], U2[N], U3[N], F1[N], F2[N], F3[N], J2[N];//守恒变量
double a[N], p[N], Ma[N];//声速，压力，马赫数
double U1_t[N], U2_t[N], U3_t[N]; //预测步
double U1_bar[N], U2_bar[N], U3_bar[N];//预测值
double F1_bar[N], F2_bar[N], F3_bar[N], J2_bar[N];
double rho_bar[N], T_bar[N],V_bar[N];
double U1_t_bar[N], U2_t_bar[N], U3_t_bar[N];//修正步
double U1_t_av[N], U2_t_av[N], U3_t_av[N];//平均值
double U1_2[N], U2_2[N], U3_2[N];//修正值
double rho_2[N], T_2[N], V_2[N];

int main()
{
	int i, j;
	ofstream file1;
	ofstream file2;
	ofstream file3;

	for (i = 0; i < N; i++)//初始条件,第0个节点到第30个节点
	{
		x[i] = i * dx;
		A[i] = 1.0 + 2.2 * (x[i] - 1.5) * (x[i] - 1.5);

		if (x[i] <= 0.5)
		{
			rho[i] = 1;
			T[i] = 1;
		}
		else if (x[i] <= 1.5)
		{
			rho[i] = 1 - 0.366 * (x[i] - 0.5);
			T[i] = 1 - 0.167 * (x[i] - 0.5);
		}
		else
		{
			rho[i] = 0.634 - 0.3879 * (x[i] - 1.5);
			T[i] = 0.833 - 0.3507 * (x[i] - 1.5);
		}
		V[i] = 0.59 / (rho[i] * A[i]);

	}

	for (j = 0; j < timestep; j++)
	{
		//初始条件,第0个节点到第30个节点
		for (i = 0; i < N; i++)
		{
			U1[i] = rho[i] * A[i];
			U2[i] = rho[i] * A[i] * V[i];
			U3[i] = rho[i] * (T[i] / (gamma - 1) + gamma / 2 * V[i] * V[i]) * A[i];

			F1[i] = U2[i];
			F2[i] = pow(U2[i], 2) / U1[i] + (gamma - 1) / gamma * (U3[i] - gamma / 2 * pow(U2[i], 2) / U1[i]);
			F3[i] = gamma * U2[i] * U3[i] / U1[i] - gamma * (gamma - 1) / 2 * pow(U2[i], 3) / pow(U1[i], 2);

		}

		//时间步长的计算		
		dt_min = CFL * dx / (sqrt(T[0]) + V[0]);
		for (i = 1; i < N; i++)
		{
			dt[i] = CFL * dx / (sqrt(T[i]) + V[i]);
			if (dt_min > dt[i])
			{
				dt_min = dt[i];
			}

		}

		//预测步
		for (i = 0; i < N - 1; i++)//计算所有内点的值，i=1~(N-2)
		{
			J2[i] = (1 / gamma) * rho[i] * T[i] * (A[i + 1] - A[i]) / dx;
			U1_t[i] = -(F1[i + 1] - F1[i]) / dx;
			U2_t[i] = -(F2[i + 1] - F2[i]) / dx + J2[i];
			U3_t[i] = -(F3[i + 1] - F3[i]) / dx;

		}

		//预测值
		for (i = 0; i < N - 1; i++)//计算所有内点的值，i=1~(N-2)
		{
			U1_bar[i] = U1[i] + U1_t[i] * dt_min;
			U2_bar[i] = U2[i] + U2_t[i] * dt_min;
			U3_bar[i] = U3[i] + U3_t[i] * dt_min;
		}

		for (i = 0; i < N - 1; i++)
		{
			rho_bar[i] = U1_bar[i] / A[i];
			V_bar[i] = U2_bar[i] / U1_bar[i];
			T_bar[i] = (gamma - 1) * (U3_bar[i] / U1_bar[i] - gamma * (U2_bar[i] / U1_bar[i]) * (U2_bar[i] / U1_bar[i]) / 2);

			F1_bar[i] = U2_bar[i];
			F2_bar[i] = pow(U2_bar[i], 2) / U1_bar[i] + (gamma - 1) / gamma * (U3_bar[i] - gamma / 2 * pow(U2_bar[i], 2) / U1_bar[i]);
			F3_bar[i] = gamma * U2_bar[i] * U3_bar[i] / U1_bar[i] - gamma * (gamma - 1) / 2 * pow(U2_bar[i], 3) / pow(U1_bar[i], 2);

		}

		//校正步
		for (i = 1; i < N - 1; i++)
		{
			J2_bar[i] = 1 / gamma * rho_bar[i] * T_bar[i] * (A[i] - A[i - 1]) / dx;
			U1_t_bar[i] = -(F1_bar[i] - F1_bar[i - 1]) / dx;
			U2_t_bar[i] = -(F2_bar[i] - F2_bar[i - 1]) / dx + J2_bar[i];
			U3_t_bar[i] = -(F3_bar[i] - F3_bar[i - 1]) / dx;
		}


		//平均值
		for (i = 1; i < N - 1; i++)
			{
				U1_t_av[i] = 0.5 * (U1_t[i] + U1_t_bar[i]);
				U2_t_av[i] = 0.5 * (U2_t[i] + U2_t_bar[i]);
				U3_t_av[i] = 0.5 * (U3_t[i] + U3_t_bar[i]);
			}

		//校正值
		for (i = 1; i < N - 1; i++)
			{
				U1_2[i] = U1[i] + U1_t_av[i] * dt_min;
				U2_2[i] = U2[i] + U2_t_av[i] * dt_min;
				U3_2[i] = U3[i] + U3_t_av[i] * dt_min;
			}

		//原始变量
		for (i = 1; i < N - 1; i++)
			{
				rho_2[i] = U1_2[i] / A[i];
				V_2[i] = U2_2[i] / U1_2[i];
				T_2[i] = (gamma - 1) * (U3_2[i] / U1_2[i] - gamma / 2 * pow(V_2[i], 2));
			}

		//边界条件
		rho_2[0] = 1;
		T_2[0] = 1;
		U1_2[0] = rho_2[0] * A[0];
		U2_2[0] = 2 * U2_2[1] - U2_2[2];
		V_2[0] = U2_2[0] / U1_2[0];
		U3_2[0] = U1_2[0] * (T_2[0] / (gamma - 1) + gamma / 2 * V_2[0] * V_2[0]);

		U1_2[N - 1] = 2 * U1_2[N - 2] - U1_2[N - 3];
		U2_2[N - 1] = 2 * U2_2[N - 2] - U2_2[N - 3];
		U3_2[N - 1] = 2 * U3_2[N - 2] - U3_2[N - 3];
		rho_2[N - 1] = U1_2[N - 1] / A[N - 1];
		T_2[N - 1] = (gamma - 1) * (U3_2[N - 1] / U1_2[N - 1] - gamma / 2 * (U2_2[N - 1] / U1_2[N - 1]) * (U2_2[N - 1] / U1_2[N - 1]));
		V_2[N - 1] = U2_2[N - 1] / U1_2[N - 1];

		for (i = 0; i < N; i++)
			{
				rho[i] = rho_2[i];
				T[i] = T_2[i];
				V[i] = V_2[i];
			}

		//压力，马赫数
		for (i = 0; i < N; i++)
			{
				p[i] = rho_2[i] * T_2[i];
				Ma[i] = V_2[i] / sqrt(T_2[i]);
			}

		if (j == 100)
		{
			file2.open("100时间步后的流场参数.txt");
			file2 << "x\t" << "rho2\t" << "V2\t" << "T2\t" << "p\t" << "Ma\t" << "U1\t" << "U2\t" << "U3\t" << "\n";
			for (i = 0; i < N; i++)
			{
				file2 << x[i] << "\t" << rho[i] << "\t" << V[i] << "\t" << T[i] << "\t" << p[i] << "\t" << Ma[i] << "\t" << U1_2[i] << "\t" << U2_2[i] << "\t" << U3_2[i] << endl;
			}
			file2.close();
		}
		if (j == 200)
		{
			file3.open("200时间步后的流场参数.txt");
			file3 << "x\t" << "rho2\t" << "V2\t" << "T2\t" << "p\t" << "Ma\t" << "U1\t" << "U2\t" << "U3\t" << "\n";
			for (i = 0; i < N; i++)
			{
				file3 << x[i] << "\t" << rho[i] << "\t" << V[i] << "\t" << T[i] << "\t" << p[i] << "\t" << Ma[i] << "\t" << U1_2[i] << "\t" << U2_2[i] << "\t" << U3_2[i] << endl;
			}
			file3.close();
		}
		}


	file1.open("1400时间步后的流场参数.txt");
	file1 << "x\t" << "rho2\t" << "V2\t" << "T2\t" << "p\t" << "Ma\t" << "U1\t" << "U2\t" << "U3\t" << "\n";
	for (i = 0; i < N; i++)
		{
			file1 << x[i] << "\t" << rho[i] << "\t" << V[i] << "\t" << T[i] << "\t" << p[i] << "\t" << Ma[i] << "\t" << U1_2[i] << "\t" << U2_2[i] << "\t" << U3_2[i] << endl;
		}
	file1.close();

	

}
