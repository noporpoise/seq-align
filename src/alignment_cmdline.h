/*
 alignment_cmdline.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_CMDLINE_HEADER_SEEN
#define ALIGNMENT_CMDLINE_HEADER_SEEN

#include "string_buffer.h"
#include "seq_file.h"

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);

void cmdline_init();
void cmdline_finish();

void cmdline_add_files(char* p1, char* p2);
int cmdline_get_num_of_file_pairs();
char* cmdline_get_file1(int i);
char* cmdline_get_file2(int i);

void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2));

#endif
