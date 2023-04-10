#include "basic_math.h"

typedef struct c_float {
	float re;
	float im;
} c_float;

inline c_float CInject(float x) {
	c_float res;
	res.re = x;
	res.im = 0;
	return res;
}

// Ok, now I get why operator overloading is nice
inline c_float CAdd(c_float x, c_float y) {
	c_float z;
	z.re = x.re + y.re;
	z.im = x.im + y.im;
	return z;
}

inline c_float CSub(c_float x, c_float y) {
	c_float z;
	z.re = x.re - y.re;
	z.im = x.im - y.im;
	return z;
}

inline c_float CMul(c_float x, c_float y) {
	c_float z;
	z.re = x.re * y.re - x.im * y.im;
	z.im = x.re * y.im + x.im * y.re;
	return z;
}

inline float CModule2(c_float z) {
	return z.re * z.re + z.im + z.im;
}

inline float CModule(c_float z) {
	return FSqrt(z.re * z.re + z.im + z.im);
}

inline c_float CConj(c_float z) {
	c_float res;
	res.re = z.re;
	res.im = - z.im;
}
 
inline c_float CDiv(c_float x, c_float y) {
	c_float z = CDiv(CMul(x, CConj(y)), CInject(CModule2(y)));
	return z;
}