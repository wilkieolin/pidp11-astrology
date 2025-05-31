#include <math.h>
#include <stdio.h>
#define PI 3.1415926

int main() {
    int num_elements = 100;
    float angles[100];
    int i = 0;
    
    if (num_elements <= 1) {
        if (num_elements == 1) {
            angles[0] = 0.0; 
        }
    } else {
        double step = (2.0 * PI) / (num_elements - 1);
        for (i = 0; i < num_elements; ++i) {
            angles[i] = (float)(-1.0 * PI + i * step);
        }
    }

    for (i = 0; i < num_elements; ++i) {
        printf("angles[%d] = %f\n", i, angles[i]);
        printf("sin(angles[%d]) = %f\n", i, sin(angles[i]));
        printf("cos(angles[%d]) = %f\n", i, cos(angles[i]));
    }

    return 0;
}