#include "basic_math.h"

struct c_float {
	float re;
	float im;

	explicit c_float(float x) : re(x), im(0) {};
	explicit c_float(float x, float y) : re(x), im(y) {};
	c_float operator+(c_float other);
	c_float operator-(c_float other);
	c_float operator*(c_float other);
	c_float operator/(c_float other);

};

inline c_float CInject(float x) {
	c_float res{ x, 0 };
	return res;
}

// Ok, now I get why operator overloading is nice
inline c_float CAdd(c_float x, c_float y) {
	c_float z{x.re + y.re, x.im + y.im};
	return z;
}

inline c_float CSub(c_float x, c_float y) {
	c_float z{ x.re - y.re, x.im - y.im };
	return z;
}

inline c_float CMul(c_float x, c_float y) {
	c_float z{ x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re };
	return z;
}

inline float CModule2(c_float z) {
	return z.re * z.re + z.im + z.im;
}

inline float CModule(c_float z) {
	return FSqrt(z.re * z.re + z.im + z.im);
}

inline c_float CConj(c_float z) {
	c_float res{ z.re, -z.im };
}
 
inline c_float CDiv(c_float x, c_float y) {
	c_float z = CDiv(CMul(x, CConj(y)), CInject(CModule2(y)));
	return z;
}

c_float c_float::operator+(c_float other) {
	return CAdd(*this, other);
}

c_float c_float::operator-(c_float other) {
	return CSub(*this, other);
}

c_float c_float::operator*(c_float other) {
	return CMul(*this, other);
}

c_float c_float::operator/(c_float other) {
	return CDiv(*this, other);
}