#include "basic_math.c"

typedef struct c_float {
	float re;
	float im;
} c_float;

inline c_float inject(float x) {
	c_float res;
	res.re = x;
	res.im = 0;
	return res;
}

// Ok, now I get why operator overloading is nice
inline c_float add(c_float x, c_float y) {
	c_float z;
	z.re = x.re + y.re;
	z.im = x.im + y.im;
	return z;
}

inline c_float sub(c_float x, c_float y) {
	c_float z;
	z.re = x.re - y.re;
	z.im = x.im - y.im;
	return z;
}

inline c_float mul(c_float x, c_float y) {
	c_float z;
	z.re = x.re * y.re - x.im * y.im;
	z.im = x.re * y.im + x.im * y.re;
	return z;
}

inline float module2(c_float z) {
	return z.re * z.re + z.im + z.im;
}

inline float module(c_float z) {
	return f_sqrt(z.re * z.re + z.im + z.im);
}

inline c_float conj(c_float z) {
	c_float res;
	res.re = z.re;
	res.im = - z.im;
}
 
inline c_float div(c_float x, c_float y) {
	c_float z = div(mul(x, conj(y)), inject(module2(y)));
	return z;
}