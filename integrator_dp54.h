#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "verification.h"
#include "integrator.h"

#define M_PIl 3.1415926535897932384626433832795028841968L

class dp54_Model : public IMathModel {
private:
	gauss* G;
	double battery;
	double Dmax;
public:

	dp54_Model(vector<double> m, gauss* G, double battery, double Dmax, vector<double> mk) {
		this->mk = mk;
		this->Dmax = Dmax;
		this->battery = battery;
		this->G = G;
		SetX0(m);
	}

	virtual double VGX() override { return G->normalRand(); }
	virtual double VGY() override { return G->normalRand(); }
	virtual double VGZ() override { return G->normalRand(); }

	vector<double> RigthParts() override {
		return X;
	}

	vector<double> f(double t, vector<double> r) override {
		const double mu = 0.012277471;
		double Ix = cos(r[8]) * sin(r[10]) * 100;
		double Iy = cos(r[6]) * sin(r[10]) * 100;
		double Iz = cos(r[6]) * cos(r[8]) * 100;
		double Ir = sin(r[6]) * sin(r[8]) * sin(r[10]) * 100;

		vector<double> rd(n);
		rd[0] = r[1];
		rd[1] = -(cos(r[6]) * sin(r[8]) * cos(r[10]) + sin(r[6]) * sin(r[10])) * U1 / mass + VGX();
		rd[2] = r[3];
		rd[3] = -(cos(r[6]) * sin(r[8]) * sin(r[10]) + sin(r[6]) * cos(r[10])) * U1 / mass + VGY();
		rd[4] = r[5];
		rd[5] = g - cos(r[6]) * cos(r[8]) * U1 / mass + VGZ() / 2;
		rd[6] = r[7];
		rd[7] = r[9] * r[11] * (Iy - Iz) / Ix - Ir / Ix * r[9] * g;
		rd[8] = r[9];
		rd[9] = r[7] * r[11] * (Iz - Ix) / Iy - Ir / Iy * r[7] * g;
		rd[10] = r[11];
		rd[11] = r[9] * r[7] * (Ix - Iy) / Iz;
		return rd;
	}
};

struct dp54_args : integratorArgs {
	double t0, hmin, hmax, h0, eps;
	double k;

public:
	dp54_args(
		IMathModel* m,
		double t0,
		double tf,
		double hmin,
		double hmax,
		double h0,
		double eps,
		double k)
		: integratorArgs(m, tf), t0(t0), hmin(hmin), hmax(hmax), h0(h0), eps(eps), k(k)
	{ }

	double getSamplingIncrement() {
		return k;
	}
};

double calcU() {
	double v = 1.0;
	double u;
	while (1.0 + v > 1.0) {
		u = v;
		v = v / 2;
	}
	return u;
}

class dp54_integrator : public integrator {
private:
	std::ofstream out; 
	const int n;
	double u;
	dp54_args args;
	double Dmax;
	double height;
public:
	dp54_integrator(dp54_args args, double d, double h) : integrator(args), args(args), n(args.model->size()) {
		out.open("results.txt");
		this->Dmax = d;
		this->height = h;
	}

	virtual bool canContinue(int t) { return t < this->tk; }

	void integration() override {
		int d = 0;
		vector<double>
			v = this->model->X0(),
			x = v,
			x0 = v;
		double
			epsMax = args.eps,
			t0 = args.t0,
			t = args.t0,
			t_out = args.t0,
			h = args.h0,
			tf = this->tk;

		vector<double> k[6];

		double
			t1 = tf,
			h_new = args.h0,
			eps;
		vector<double> x1(n);

		int i;
		int j = 0;
		double epsD = 0.9;
		vector<double> temp(12); //текущий
		vector<double> prev(12);
		for (i = 0; i < n; i++) {
			temp[i] = x[i];
		}
		double r = D(v, Dmax);
		while (d < dMax)
		{
			if ((temp[4] <= (height - height / 5)) && (abs(temp[0] - args.model->mk_()[0]) >= 0.1)) {
				h = h_new;
				if (t + h >= tf)
					h = tf - t;
				oneStep(h, t, x, epsMax, t1, x1, k, eps); temp = x;
				h_new = h / (max(0.1, min(5.0, pow(eps / epsMax, 0.2) / 0.9)));

				if (eps > epsMax)
					continue;

				while ((t_out < t + h) && (t_out <= t1))
				{
					double
						theta = (t_out - t) / h,
						b[6];

					double T = theta;
					b[0] = T * (1 + T * (-1337 / 480.0 + T * (1039 / 360.0 + T * (-1163 / 1152.0))));
					b[1] = 0;
					b[2] = 100.0 * T * T * (1054 / 9275.0 + T * (-4682 / 27825.0 + T * (379 / 5565.0))) / 3.0;
					b[3] = -5.0 * T * T * (27 / 40.0 + T * (-9 / 5.0 + T * (83 / 96.0))) / 2.0;
					b[4] = 18225.0 * T * T * (-3 / 250.0 + T * (22 / 375.0 + T * (-37 / 600.0))) / 848.0;
					b[5] = -22.0 * T * T * (-3 / 10.0 + T * (29 / 30.0 + T * (-17 / 24.0))) / 7.0;

					vector<double> xOut(n);
					for (int mi = 0; mi < n; mi++) {
						xOut[mi] = x[mi];
						for (int mj = 0; mj < 6; mj++) {
							xOut[mi] += b[mj] * k[mj][mi];
						}
					}

					SetX0(xOut);
					r = D(xOut, Dmax);
					for (i = 0; i < n; i++) {
						temp[i] = xOut[i];
					}
					cout << xOut[0] << "        " << xOut[2] << "                  " << xOut[4] << endl;
					j++;
					out << xOut[0] << "   " << xOut[2] << "   " << xOut[4] << endl;
					t_out += this->args.getSamplingIncrement();
				}
				for (i = 0; i < n; i++) {
					x[i] = x1[i];
				}
				t += h;
				
			}

			if ((temp[4] >= (height - height / 5) && (temp[4] < height))) {
				temp[4] = temp[4] + 0.1 * height;
				temp[5] = temp[5] - 0.4;
				if (args.model->mk_()[0] > temp[0]) {
					temp[0] = temp[0] + 0.2;
				}
				else temp[0] = temp[0] - 0.2;
				for (i = 0; i < n; i++) {
					prev[i] = temp[i];
				}
				SetX0(temp);
				out << temp[0] << "   " << temp[2] << "   " << temp[4] << endl;
				r = D(temp, Dmax);
			}
			if ((temp[4] >= height)) {
				if ((temp[0] >= (args.model->mk_()[0] - 1)) && (temp[0] <= (args.model->mk_()[0] + 1))) {
					while (abs(temp[4] - args.model->mk_()[2]) > 1) {
						if (args.model->mk_()[0] > temp[0]) {
							temp[0] = temp[0] + 0.1;
						}
						else temp[0] = temp[0] - 0.1;
						if (args.model->mk_()[1] > temp[2]) {
							temp[2] = temp[2] + 0.1;
						}
						else temp[2] = temp[2] - 0.1;

						if (args.model->mk_()[2] > temp[4]) {
							temp[4] = temp[4] + 2;
						}
						else temp[4] = temp[4] - 2;

						SetX0(temp);
						cout << temp[0] << "        " << temp[2] << "                  " << temp[4] << endl;
						out << temp[0] << "   " << temp[2] << "   " << temp[4] << endl;
						r = D(temp, Dmax);
					}
				}
				else
				{
					if (args.model->mk_()[0] > temp[0]) {
						temp[0] = temp[0] + 0.2;
					}
					else temp[0] = temp[0] - 0.2;
					temp[2] = (temp[0] * prev[2] - temp[0] * args.model->mk_()[1] - args.model->mk_()[0] * prev[2] + args.model->mk_()[2] * prev[0]) / (prev[0] - args.model->mk_()[0]);
					SetX0(temp);
					r = D(temp, Dmax);
					cout << temp[0] << "        " << temp[2] << "                  " << temp[4] << endl;
					out << temp[0] << "   " << temp[2] << "   " << temp[4] << endl;
				}

			}
			d++;
		}
	}


	double getSamplingIncrement() {
		return 1;
	}

private:
	void oneStep(
		double h,
		double t0,
		vector<double> r0,
		double epsMax,
		double& t1,
		vector<double>& r1,
		vector<double> k[6],
		double& erra);

	void f(double t, vector<double> r, vector<double>& rd);
};

void dp54_integrator::f(double t, vector<double> x0, vector<double>& x1) {
	x1 = model->f(t, x0);
}

double _max(double d1, double d2, double d3, double d4) {
	return max(max(d1, d2), max(d3, d4));
}

void dp54_integrator::oneStep(
	double h,
	double t0,
	vector<double> x0,
	double epsMax,
	double& t1,
	vector<double>& x1,
	vector<double> k[6],
	double& eps)
{
	vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n);
	double
		a2 = 1.0L / 5,
		b21 = 1.0L / 5;
	double
		a3 = 3.0L / 10,
		b31 = 3.0L / 40,
		b32 = 9.0L / 40;
	double
		a4 = 4.0L / 5,
		b41 = 44.0L / 45,
		b42 = -56.0L / 15,
		b43 = 32.0L / 9;
	double
		a5 = 8.0L / 9,
		b51 = 19372.0L / 6561,
		b52 = -25360.0L / 2187,
		b53 = 64448.0L / 6561,
		b54 = -212.0L / 729;
	double
		a6 = 1.0L,
		b61 = 9017.0L / 3168,
		b62 = -355.0L / 33,
		b63 = 46732.0L / 5247,
		b64 = 49.0L / 176,
		b65 = -5103.0L / 18656;
	double
		a7 = 1.0L,
		b71 = 35.0L / 384,
		b73 = 500.0L / 1113,
		b74 = 125.0L / 192,
		b75 = -2187.0L / 6784,
		b76 = 11.0L / 84;
	double
		c1 = 35.0L / 384,
		c3 = 500.0L / 1113,
		c4 = 125.0L / 192,
		c5 = -2187.0L / 6784,
		c6 = 11.0L / 84;
	double tk0, tk1, tk2, tk3, tk4, tk5, tk6;
	vector<double> yk0(n), yk1(n), yk2(n), yk3(n), yk4(n), yk5(n), yk6(n);
	vector<double> rz(n);
	double s;
	int i;
	tk0 = t0;
	for (i = 0; i < n; i++)
		yk0[i] = x0[i];
	f(tk0, yk0, k1);
	for (i = 0; i < n; i++)
		k1[i] *= h;
	tk1 = t0 + a2 * h;
	for (i = 0; i < n; i++)
		yk1[i] = x0[i] + b21 * k1[i];
	f(tk1, yk1, k2);
	for (i = 0; i < n; i++)
		k2[i] *= h;
	tk2 = t0 + a3 * h;
	for (i = 0; i < n; i++)
		yk2[i] = x0[i] + b31 * k1[i] + b32 * k2[i];
	f(tk2, yk2, k3);
	for (i = 0; i < n; i++)
		k3[i] *= h;
	tk3 = t0 + a4 * h;
	for (i = 0; i < n; i++)
		yk3[i] = x0[i] + b41 * k1[i] + b42 * k2[i] + b43 * k3[i];
	f(tk3, yk3, k4);
	for (i = 0; i < n; i++)
		k4[i] *= h;
	tk4 = t0 + a5 * h;
	for (i = 0; i < n; i++)
		yk4[i] = x0[i] + b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i];
	f(tk4, yk4, k5);
	for (i = 0; i < n; i++)
		k5[i] *= h;
	tk5 = t0 + a6 * h;
	for (i = 0; i < n; i++)
		yk5[i] = x0[i] + b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i];
	f(tk5, yk5, k6);
	for (i = 0; i < n; i++)
		k6[i] *= h;
	tk6 = t0 + a7 * h;
	for (i = 0; i < n; i++)
		yk6[i] = x0[i] + b71 * k1[i] + b73 * k3[i] + b74 * k4[i] + b75 * k5[i] + b76 * k6[i];
	f(tk6, yk6, k7);
	for (i = 0; i < n; i++)
		k7[i] *= h;
	for (i = 0; i < n; i++)
		x1[i] = x0[i] + c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i];
	// e calculated as difference between x1^ and x1
	double
		e1 = 71.0L / 57600,
		e3 = -71.0L / 16695,
		e4 = 71.0L / 1920,
		e5 = -17253.0L / 339200,
		e6 = 22.0L / 525,
		e7 = -1.0L / 40;
	vector<double> e = { e1, 0, e3, e4, e5, e6, e7 };
	k[0] = k1;
	k[1] = k2;
	k[2] = k3;
	k[3] = k4;
	k[4] = k5;
	k[5] = k6;
	{
		// calculate erra
		eps = 0.0;
		for (i = 0; i < n; i++)
		{
			double d = 0;
			for (int j = 0; j < 6; j++) {
				d += e[j] * k[j][i];
			}
			eps += pow((h * d) / _max(pow(10.0, -5.0), abs(x1[i]), abs(x0[i]), 2 * u / epsMax), 2.0);
		}
		eps = sqrt(eps / n);
	}
}