#include "func.h"

int strToInt(char* s) {
    int res = 0;
    while (*s != '\0') {
        res = res * 10 + *s - '0';
        s++;
    }
    return res;
}

