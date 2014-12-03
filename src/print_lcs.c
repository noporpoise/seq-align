#include "print_lcs.h"

void print_lcs(alignment_t *result)
{
  int i = 0;
  printf("LCS: \n");
  for(; i < strlen(result->result_a); i++)
    if(result->result_a[i] == result->result_b[i])
      printf("%c", result->result_a[i]);
  printf("\n");
}