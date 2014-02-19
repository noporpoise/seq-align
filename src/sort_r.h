/* Isaac Turner 23 Jan 2014 Public Domain */
#include <stdlib.h>

/*
sort_r function to be exported.

Parameters:
  base is the array to be sorted
  nel is the number of elements in the array
  width is the size in bytes of each element of the array
  compar is the comparison function
  arg is a pointer to be passed to the comparison function

void sort_r(void *base, size_t nel, size_t width,
            int (*compar)(const void *a1, const void *a2, void *aarg),
            void *arg);
*/

#if defined(__MINGW32__) || defined(__OpenBSD__) || defined(AMIGA) || \
    defined(__gnu_hurd__) || (__GLIBC__ == 2 && __GLIBC_MINOR__ < 8)
  #define QSORT_WITH_NESTED_FUNCTIONS 1
#endif

#ifdef QSORT_WITH_NESTED_FUNCTIONS

static inline void sort_r(void *base, size_t nel, size_t width,
                          int (*compar)(const void *a1, const void *a2, void *aarg),
                          void *arg)
{
  int nested_cmp(const void *a, const void *b)
  {
    return compar(a, b, arg);
  }

  qsort(base, nel, width, nested_cmp);
}

#else

#if (defined _GNU_SOURCE || defined __GNU__ || defined __linux__)

typedef int(* __compar_d_fn_t)(const void *, const void *, void *);
extern void qsort_r (void *__base, size_t __nmemb, size_t __size,
                     __compar_d_fn_t __compar, void *__arg)
  __nonnull ((1, 4));

#else /* Not GNU linux */

struct sort_r_data
{
  void *arg;
  int (*compar)(const void *a1, const void *a2, void *aarg);
};

static inline int sort_r_arg_swap(void *s, const void *aa, const void *bb)
{
  struct sort_r_data *ss = (struct sort_r_data*)s;
  return (ss->compar)(aa, bb, ss->arg);
}

#endif

static inline void sort_r(void *base, size_t nel, size_t width,
                          int (*compar)(const void *a1, const void *a2, void *aarg),
                          void *arg)
{
  #if (defined _GNU_SOURCE || defined __GNU__ || defined __linux__)

    qsort_r(base, nel, width, compar, arg);

  #elif (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || \
         defined __FREEBSD__ || defined __BSD__ || \
         defined OpenBSD3_1 || defined OpenBSD3_9)

    struct sort_r_data tmp;
    tmp.arg = arg;
    tmp.compar = compar;
    qsort_r(base, nel, width, &tmp, &sort_r_arg_swap);

  #elif (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)

    struct sort_r_data tmp;
    tmp.arg = arg;
    tmp.compar = compar;
    qsort_s(base, nel, width, &sort_r_arg_swap, &tmp);

  #else
    #error Cannot detect operating system
  #endif
}

#endif /* !QSORT_WITH_NESTED_FUNCTIONS */
