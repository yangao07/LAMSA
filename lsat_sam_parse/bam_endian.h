#ifndef BAM_ENDIAN_H
#define BAM_ENDIAN_H

#include <stdint.h>

// XXX s
static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}
// XXX e


#endif
