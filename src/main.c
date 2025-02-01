#include "stdio.h"
#include "vinc.h"

int main() {
    // Testing with one value
    struct ResultTrans res = trans(40.99698, 46.0, 9.20127, 10.0);
    printf("%f - %f\n", res.x, res.y);
    return 0;
}