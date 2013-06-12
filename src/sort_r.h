
#include <stdlib.h>

/* sort_r function to be exported */
void sort_r(void *base, size_t nel, size_t width,
            int (*compar)(const void *a1, const void *a2, void *aarg), void *arg);
