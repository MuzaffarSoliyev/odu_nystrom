#include "task.h"
#include "funkr.h"

double func(double x,double y1,double y2){
	return f(x) - q(x)*y1 - p(x)*y2;
}

void Runge_Kutta(double *y1,double *y2, double y0, double x0, double h) {
	double k1, k2, k3, k4, x;
	double m1, m2, m3, m4;
	for (int i = 0; i < 4; i++){
		x = x0 + h*i;

		m1 = h*y2[i];
		k1 = h*func(x, y1[i], y2[i]);

		m2 = h*(y2[i] + k1 / 2.0);
		k2 = h*func(x + h / 2.0, y1[i] + m1 / 2.0, y2[i] + k1 / 2.0);

		m3 = h*(y2[i] + k2 / 2.0);
		k3 = h*func(x + h / 2.0, y1[i] + m2 / 2.0, y2[i] + (k2 / 2.0));

		m4 = h*(y2[i] + k3);
		k4 = h*func(x + h, y1[i] + m3, y2[i] + k3);

		y1[i + 1] = y1[i] + (m1 + 2 * (m2 + m3) + m4) / 6.0;
		y2[i + 1] = y2[i] + (k1 + 2 * (k2 + k3) + k4) / 6.0;
	}
	
	
}

double Nystrom(double *y1, double *y2, double y0, double x0, double h, int n,double f1,double val) {
	double k1_y1,k1_y2,k2_y1,k2_y2,k3_y1,k3_y2, k4_y2, k4_y1;
	double delta_1, delta_2;
	double temp1, temp2,temp3,temp4,temp5,temp6;

	y1[0] = y0;
	y2[0] = val;

	Runge_Kutta(y1,y2,y0, x0, h); //–унге- утта 4 пор€дка

	for (int i = 3; i < n; i++) {
        /* явный Ќюстрем. ѕантелеев ќбыкновенные дифференциальные уравнени€ в примерах и задачах. стр 350*/
		/* ’айрер. стр 329 */

		k1_y1 = y2[i];
		k1_y2 = func(x0 + i*h, y1[i], y2[i]);
		
		k2_y1 = y2[i-1];
		k2_y2 = func(x0 + (i-1)*h, y1[i-1], y2[i-1]);

		k3_y1 = y2[i-2];
		k3_y2 = func(x0 + (i-2)*h, y1[i-2], y2[i-2]);

		k4_y1 = y2[i-3];
		k4_y2 = func(x0 + (i-3)*h, y1[i-3], y2[i-3]);
		
		delta_1 = h*(8 * k1_y1 - 5 * k2_y1 + 4 * k3_y1 - k4_y1) / 3;
		delta_2 = h*(8 * k1_y2 - 5 * k2_y2 + 4 * k3_y2 - k4_y2) / 3;
		
		y1[i+1] = y1[i-1] + delta_1;
		y2[i+1] = y2[i-1] + delta_2;
	}

	return y1[n] - f1;
}

int iteration(double a, double b, double f0, double f1,double *y_1, double *y_2, int n,double h){
	double test_x_1 = 0;
	double test_x_2 = 0;
	double test_x_3 = 0;

	double func_2 = 0;
	double func_1 = 0;
	double func_3 = 0;

	test_x_1 = -10;
	test_x_2 = 10;

	double border = 0.00000001;//подгон€етс€ под задачу

	while (1) {
		func_1 = Nystrom(y_1, y_2, f0, a, h, n, f1, test_x_1);
		func_2 = Nystrom(y_1, y_2, f0, a, h, n, f1, test_x_2);

		if (func_1*func_2 > 0) {
			break;
		}

		test_x_3 = (test_x_1 + test_x_2) / 2;
		func_3 = Nystrom(y_1, y_2, f0, a, h, n, f1, test_x_3);

		if (func_3*func_1 > 0) {
			test_x_1 = test_x_3;
		}

		if (func_3*func_2 > 0) {
			test_x_2 = test_x_3;
		}

		printf("%lf %lf %lf\n", func_1, func_2, func_3);

		if ((fabs(func_1) < border) || (fabs(func_2) < border)) {
			break;
		}
	}
}

int odu(double a, double b, double f0, double f1, double e, double *y_1,double *y_2, int n) {
	double *test_1 = (double*)odu_memsize(1000);
	double *test_2 = (double*)odu_memsize(1000);

	int max_iter = 0;

	double h1 = (b - a) / n;
	double h2 = (b - a) / (2 * n);

	int n1 = n;
	int n2 = 2 * n1;
	int index1 = 0;
	int index2 = 0;
	int flag = 1;

	while (max_iter < 10){
		iteration(a, b, f0, f1, y_1, y_2, n1, h1);
		for (int i = 0; i < n1; i++){
			if (i%(n1/n) == 0) {
				test_1[index1] = y_1[i];
				index1 = index1 + 1;
			}
		}

		memset(y_1, 0, sizeof(y_1));
		memset(y_2, 0, sizeof(y_2));

		iteration(a, b, f0, f1, y_1, y_2, n2, h2);
		for (int i = 0; i < n2; i++) {
			if (i%(n2/n) == 0) {
				test_2[index2] = y_1[i];
				index2 = index2 + 1;
			}
		}

		for (int i = 0; i < n; i++) {
			if (fabs(test_1[i] - test_2[i]) > e){
				flag = 0;
				break;
			}
		}

		if (flag==1){
			printf("Convergent!");
			memset(y_1, 0, sizeof(y_1));
			for (int i = 0; i < n; i++) {
				y_1[i] = test_1[i];
			}
			getchar();
			return n;
		}

		memset(y_1, 0, sizeof(y_1));
		memset(y_2, 0, sizeof(y_2));
		memset(test_1, 0, sizeof(test_1));
		memset(test_2, 0, sizeof(test_2));

		n1 = n2;
		n2 = 2 * n2;
		h1 = h2;
		h2 = h2 / 2;
		max_iter++;
		index1 = 0;
		index2 = 0;
		flag = 1;
	}

	return n1;
}

