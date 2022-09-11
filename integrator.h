#pragma once
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <fstream>
#include <random>
#include <conio.h>
#include "UI.h"

using namespace std;
const double g = 9.8;
const int m = 2;
const int a = 1;
const int n = 12;
const double Cx = 0.47;
const int dMax = 10000;
const double S = 0.00785;
const double mass = 1;
const double U1 = 1;

class gauss {
private:
	double mean;
	double stddev;

	double uniformRand()
	{
		return static_cast<double>(rand()) / RAND_MAX;
	}

public:

	double normalRand()
	{
		return sqrt(-2.0 * log(uniformRand())) * cos(2.0 * 3.14 * uniformRand()) * stddev + mean;
	}

	gauss(double mean, double stddev) {
		this->mean = mean;
		this->stddev = stddev;
	}
};

class IMathModel {
protected:
	double Dmax;
	vector<double> mk;
	vector<double> X;
	vector<vector<double>> Vectors;
public:
	virtual vector<double> RigthParts() = 0;

	virtual void SetX0(vector<double> x0) {
		X = x0;
		Vectors.push_back(x0);
	}

	virtual double VGX() { return 0; }
	virtual double VGY() { return 0; }
	virtual double VGZ() { return 0; }

	vector<double> X0() { return X; }
	vector<double> mk_() { return mk; }
	vector<vector<double>> X_all() { return Vectors; }

	int size() { return X.size(); }

	void print() {
		int n = X.size();
		for (int i = 0; i < n; i++) {
			cout << X[i] << " ";
		}
		cout << endl;
	};

	virtual vector<double> f(double t, vector<double> r) { return r; }
};


struct integratorArgs {
	integratorArgs(IMathModel* m, double tk) {
		model = m;
		this->tk = tk;
	}
	IMathModel* model = nullptr;
	double tk;
};

class integrator {
protected:
	IMathModel* model;
	double tk;
public:

	integrator(integratorArgs args) : tk(args.tk), model(args.model) { }

	~integrator() {
		cout << "class is destroyed" << endl;
	}

	vector<double> RigthParts() {
		return model->RigthParts();
	}

	void SetX0(vector<double> x0) {
		model->SetX0(x0);
	}

	void print() {
		model->print();
	}

	virtual void integration() = 0;
};

