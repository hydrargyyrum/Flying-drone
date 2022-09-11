#define _CRT_SECURE_NO_WARNINGS

#include "integrator.h"
#include "integrator_dp54.h"
#include "verification.h"
#include <functional>

IMathModel* test_dp54() {
	double
		t0 = 0,
		tf = 5000,
		h0 = 0.001,
		epsMax = 0.001,
		k = 0.01,
		battery = 100, Dmax = 5000;
	gauss G(-1, 3);
	vector<double> mk = { -10,10,1 };
	//green - x; red- y, blue - z
	vector<double> m = { 1, 2 , 0,1,0,0,0,-10,0,0,-10,0 };
	IMathModel* modelA = new dp54_Model(m, &G, battery, Dmax, mk);
	dp54_args args(modelA, t0, tf, 0, 0, h0, epsMax, k);
	dp54_integrator c1(args,Dmax,Dmax/50);
	if (checkBattery(Dmax, battery, mk)&&(safetyDistanceForFligth(m,mk))) {
		c1.integration();
	}
	else cout << "I don't have enough charge for the fligth or unsafe fligth distance set" << endl;
	return modelA;
}

void attractor_show(IMathModel* model, function<vector<double>(vector<double>)> map) {
	attractor_show_args args_attr;
	args_attr.example_to_show = false;
	args_attr.example_continiously_draw = false;
	args_attr.animate = true;
	auto v = model->X_all();
	cout << "size : " << v.size() << endl;
	vector<vector<double>> points(v.size());
	std::transform(v.begin(), v.end(), points.begin(), map);
	args_attr.points = points;
	args_attr.print_dots = false;
	attractor_show(args_attr);
}

void main() {
	srand(time(nullptr));
	attractor_show(test_dp54(), [](vector<double> v) -> vector<double> {
		return vector<double>{ v[0] *4, v[2]*4 , v[4]*0.5};
		});
	std::cin.get();
}