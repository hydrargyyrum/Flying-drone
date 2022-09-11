#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
using namespace std;

static bool checkBattery(double Dmax, double battery, vector<double> v) {
	if ((sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < Dmax) && 
		(sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * 0.015 < battery)) return true;
	else return false;
}

static double D(vector<double> v, double D) {
	return abs(sqrt(v[0] * v[0] + v[2] * v[2] + v[4] * v[4]) - D);
}

static bool safetyDistanceForFligth(vector<double> v1, vector<double> v2) {
	if (sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[2] - v2[1]) * 
		(v1[2] - v2[1]) + (v1[4] - v2[2]) * (v1[4] - v2[2])) < 10) return false;
	else return true;
}