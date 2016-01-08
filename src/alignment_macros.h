/*
 alignment_macros.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_MACROS_HEADER_SEEN
#define ALIGNMENT_MACROS_HEADER_SEEN

#define ARR_2D_INDEX(width,i,j) (((unsigned long)(j)*(width)) + (i))
#define ARR_LOOKUP(arr,width,i,j) arr[ARR_2D_INDEX((width),(i),(j))]
#define ARR_2D_X(arr_index, arr_width) ((arr_index) % (arr_width))
#define ARR_2D_Y(arr_index, arr_width) ((arr_index) / (arr_width))

#define QUOTE(str) #str

#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))
#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))
#define MAX4(w,x,y,z) MAX2(MAX3(w, x, y), z)

#define ABSDIFF(a,b) ((a) > (b) ? (a)-(b) : (b)-(a))

#endif /* ALIGNMENT_MACROS_HEADER_SEEN */
