#include <mpi.h>
#include <stdio.h>
#include <windows.h>

double f(double x)
{
	return 2 * x * (x + 0.2) + 0.4;
}
double psi1(double t)
{
	return 2 * t + 0.4;
}
double psi2(double t)
{
	return 1.36;
}

int main(int argc, char* argv[])
{
	system("chcp 1251>nul");
	int size, rank, msgtag = 14; MPI_Status status;
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) { return 1; }
	if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) { MPI_Finalize(); return 2; }
	if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) { MPI_Finalize(); return 3; }
	//====================================================================================

	//Simulation parameters
	double t = 0.0, tmax = 100.0, tau = 0.01;
	double x = 0.0, Len = 300.0, h = 0.001;
	double a = 1.0, R = (a * tau) / (h * h);
	double tn, tk, delta;
	double tmp;

	//Data arrays
	int N = (int)(Len / h) + 1;
	double* U = new double[N];
	double* L = new double[N];
	double* M = new double[N];

	//Initial state of rod
	U[0] = psi1(t);
	x = h;
	for (int i = 1; i < N - 1; i++) U[i] = f(x);
	U[N - 1] = psi2(t);

	//Initial values of L & M coefficients
	L[0] = 0;
	M[0] = psi1(t);

	//Calc L coefficient values
	for (int i = 1; i < N; i++) L[i] = R / (1 + 2 * R - R * L[i - 1]);

	//Start timer
	tn = MPI_Wtime();

	//Main calculation loop
	do
	{
		//Straight step of sweep method
		U[0] = psi1(t + tau);
		for (int i = 1; i < N; i++) M[i] = (U[i] + R * M[i - 1]) / (1 + 2 * R - R * L[i - 1]);
		U[N - 1] = psi2(t + tau);

		//Reverse step of sweep method
		for (int i = N - 2; i > 0; i--) U[i] = L[i] * U[i + 1] + M[i];

		//Time step
		t += tau;
		printf_s("Time: %f\n",t);
	} while (t <= tmax);

	//Stop timer
	tk = MPI_Wtime();
	delta = tk - tn;

	//Output
	FILE* f;
	fopen_s(&f, "SEQ_Result.txt", "w");
	if (f)
	{
		fprintf_s(f, "Time: %4.4f\n", delta);
		for (int i = 0; i < N; i++) fprintf_s(f, "U[%d] = %f\n", i, U[i]);
		fclose(f);
		printf_s("Output is in the file.\n");
	}
	else printf_s("File error!\n");

	//Free memory
	delete[]U;
	delete[]L;
	delete[]M;

	//====================================================================================
	MPI_Finalize();
	system("pause>nul");
	return 0;
}