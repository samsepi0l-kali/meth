#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

int main() {
    const int n = 3;
    const int p = 3;

    // Variant 1 matrix
    double A[n][n] = {
        { 3, -7,  -1},
        { -9,  -8, -7},
        { 5, 2,  2}
    };

    // Variant 1 vector
    double b[n] = {1, 5, -3};

    double x[p];   // Vector p-norms
    double y;      // Vector infinity norm
    double u;      // Matrix 1-norm
    double v;      // Matrix 2-norm
    double w;      // Matrix infinity norm

    // Calculate vector norms for p = 1, 2, 3
    for (int j = 1; j <= p; j++) {
        x[j - 1] = 0;

        for (int i = 0; i < n; i++) {
            double h = abs(b[i]);
            x[j - 1] += pow(h, j);
        }

        x[j - 1] = pow(x[j - 1], 1.0 / j);
    }

    // Vector infinity norm
    y = 0;

    for (int i = 0; i < n; i++) {
        y = max(y, abs(b[i]));
    }

    // Matrix 1-norm:
    // maximum absolute column sum
    u = 0;

    for (int j = 0; j < n; j++) {
        double h = 0;

        for (int i = 0; i < n; i++) {
            h += abs(A[i][j]);
        }

        u = max(u, h);
    }

    // Matrix 2-norm according to the lecture program:
    // square root of the sum of squares
    v = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            v += A[i][j] * A[i][j];
        }
    }

    v = sqrt(v);

    // Matrix infinity norm:
    // maximum absolute row sum
    w = 0;

    for (int i = 0; i < n; i++) {
        double h = 0;

        for (int j = 0; j < n; j++) {
            h += abs(A[i][j]);
        }

        w = max(w, h);
    }

    // Output
    cout << fixed << setprecision(6);

    cout << "Vector norms:" << endl;
    cout << "||b||1   = " << x[0] << endl;
    cout << "||b||2   = " << x[1] << endl;
    cout << "||b||3   = " << x[2] << endl;
    cout << "||b||inf = " << y << endl;

    cout << endl;

    cout << "Matrix norms:" << endl;
    cout << "||A||1   = " << u << endl;
    cout << "||A||2   = " << v << endl;
    cout << "||A||inf = " << w << endl;

    return 0;
}