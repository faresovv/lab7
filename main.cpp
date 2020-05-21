#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void adams_method(double h, double *y, double *dy) {
    int n = (2+h)/h;

    double x[n];
    for(int i = 0; i < n; i++)
        x[i] = i*h;

    y[0] = 1;
    dy[0] = 0;
    double d2y_p = 16;
    y[1] = y[0] + h*dy[0] + h*h*d2y_p/2.;
    dy[1] = dy[0] + h*d2y_p + h*h*(4-16*d2y_p)/2.;
    double d2y_c = sin(4*x[1])*exp(x[1]) - 16*dy[1] + 16*y[1];
    for(int i = 2; i < n; i++) {
        y[i] = y[i-1] + h/2.*(3*dy[i-1] - dy[i-2]);
        dy[i] = dy[i-1] + h/2.*(3*d2y_c - d2y_p);
        d2y_p = sin(4*x[i-1])*exp(x[i-1]) - 16*dy[i-1] + 16*y[i-1];
        d2y_c = sin(4*x[i])*exp(x[i]) - 16*dy[i] + 16*y[i];
    }
}

double pravilo_runge(double eps = 0.01) {
    double h = 1;
    bool conv = true;
    do {
        int n = (2+h)/h;
        double y[n], dy[n];
        adams_method(h, y, dy);
        h /= 2;
        double y1[2*n], dy1[2*n];
        adams_method(h, y1, dy1);

        double delta_y = abs(y[0] - y1[0]);
        double delta_dy = abs(dy[0] - dy1[0]);

        for(int i = 0; i < n; i++) {
            if(delta_y < abs(y[i] - y1[2*i])){
                delta_y = abs(y[i] - y1[2*i]);
            }
            if(delta_dy < abs(dy[i] - dy1[2*i])) {
                delta_dy = abs(dy[i] - dy1[2*i]);
            }
        }
        conv = delta_y >= eps || delta_dy >= eps;
    } while(conv);
    return h;
}

int main()
{
    double h = pravilo_runge();
    int n = (2+h)/h;
    double y[n], dy[n];
    adams_method(h, y, dy);

    ofstream out("ans.dat", ios_base::out);

    for(int i = 0; i < n; i++){
        out << y[i] << "," << dy[i] << "," << i*h << endl;
    }
    out.close();
    return 0;
}
