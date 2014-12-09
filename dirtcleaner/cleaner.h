#ifndef _CLEANER_H_
#define _CLEANER_H_

#include "stdafx.h"
#include <stdio.h>
#include <assert.h>
#include <vd2/plugin/vdvideofilt.h>
#include <vector>

void log(const char* str);
void logn(const char* str, int n);
#define SHOW(x) logn(#x, x);

template<int W, typename T = BYTE>
class __declspec( align(16)) MonoBlock
{
public:
	T data[W*W];

	/*void split(MonoBlock<W/2, T> **subs);	
	void split(MonoBlock<W/2, T> *subs);	
	void combine(MonoBlock<W/2, T> **subs);
	void combine_f(MonoBlock<W/2, float> **subs); // T must be BYTE
	void print(const char *name);
	void print_f(const char *name);

	void halve(MonoBlock<W/2, T> &res);
	void merge_f(MonoBlock<W,float> &a, MonoBlock<W,float> &b, float k); //T must be float*/
	MonoBlock() { /*memset(data, 255, W*W);*/ } //for debug
	//MonoBlock(const MonoBlock<W, BYTE> &a); //copy
	//void round(MonoBlock<W,float> &a);
};

template<class T>
class TVec
{
public:
	T x, y;

	TVec(T vx, T vy) : x(vx), y(vy) {}
	TVec() : x(0), y(0) {}	
	TVec operator+(const TVec &a) const	{	return TVec(x + a.x, y + a.y);	}
	TVec operator-(const TVec &a) const	{	return TVec(x - a.x, y - a.y);	}
	TVec& operator+=(const TVec &a)		{	x += a.x; y += a.y;	return *this;	}	
	TVec plus(const TVec &a)	const	{	return TVec(x + a.x, y + a.y); }
	TVec& move(const TVec &a)			{	x += a.x; y += a.y;	return *this;	}
	//void print(const char* name)	const	{	printf("%s(%d, %d) ", name, x, y);   }
	bool operator==(const TVec &a) const { return x==a.x && y==a.y;  }
	TVec<float> operator*(float k) const { return TVec<float>(x*k, y*k); }
	TVec div2() const { return TVec(x / 2, y / 2); }
};

typedef TVec<int> Vec;
typedef TVec<float> FVec;

class Plane
{
public: 
	Plane() : stride(0), width(0), height(0), border(32), data(NULL), offset(0), ownData(false) {}
	~Plane();
	void create(int w, int h);
	void destroy();
	void copyFrom(BYTE* data, int w, int h, int pitch);
	//void alias(Plane &other); //be an alias of some other plane

	//template<int W>	void writeBlock(int x0, int y0, MonoBlock<W> &block);
	//template<int W>	void readBlock(int x0, int y0, MonoBlock<W> &block);
	template<int W> void readBlock_f(Vec v, MonoBlock<W, float> &block);

	/*void halveRect(Vec<VP_HIRES> v, int bw, int bh, int bstride, BYTE *buf);
	template<int W>	PBlock<W> makePBlock() { return PBlock<W>(stride); }
	template<int W>	void getPBlock(Vec<VP_LOWRES> v, PBlock<W> &block);
	template<int W> void halveBlockHP(Vec<VP_HALF> v, MonoBlock<W,float> &block);
	template<int W> void halveBlockQP(Vec<VP_QUARTER> v, MonoBlock<W,float> &block);
	template<int W> void readBlockQP(Vec<VP_QUARTER> v, MonoBlock<W,float> &block);*/
	void swap(Plane &other);
	void makeBorder(int startRow, int nRows);
	//__int64 diff(Plane &other);
	//void mkHalfLuma(Plane &halfLuma, int startRow, int nRows);

	//template<int W> void halveBlockVertQP(int x, int y, int dy, const F32vec4 * vsrc, F32vec4 *vtmp, bool oneMoreCol); //aligned float tmp[W * (W*2+4 + 1)];
	//template<int W> void halveBlockVertQP2(int x, int y, F32vec4 *vtmp, bool oneMoreCol); //aligned float tmp[W * (W*2+4 + 1)];

	//template<int W> void halveBlockVertHP(int x, int y, const F32vec4 * __restrict vsrc, float * __restrict tmp);
	//template<int W> void readRect_f(int x, int y, float * __restrict src); //[src_stride * (src_stride + 1)]  where src_stride = W*2+4

	BYTE* pixelPtr(int y, int x);
	//void get16floats(int x, int y, F32vec4 *v); //v[at least 3]

	//friend class RgbPlane;
	BYTE *data; //16 bytes aligned, includes border
	int height;
//protected: //for gym
	int width;
	int stride, border, offset;
	bool ownData;
};

class YV12Plane {
public:
	Plane Y,U,V;
	void create(int w, int h);
	void destroy();
	~YV12Plane();
	void swap(YV12Plane &other);
	void copyFrom(const VDXPixmap &src);
};

typedef std::vector< std::vector< Vec > > VecMatrix;

class Cleaner {
public:
	Cleaner() : fn(0), nbx(0), nby(0), inited(0) {}

	void init(int w, int h);
	void process(const VDXPixmap &src, const VDXPixmap &dst, int nFrame);

	int inited;

protected:
	//Vec motionSearch(int bx, int by, MonoBlock<8, float> &srcBlock, Plane &plane); // => offset vector
	Vec getMVp(int bx, int by);
	Vec getMVpCenter(int bx, int by);
	Vec getMVn(int bx, int by);
	Vec getMVnCenter(int bx, int by);
	Vec getMVCenter(int bx, int by, bool prev);
	void flowBlock(int bx, int by, bool prev, BYTE* yv12block); // yv12block [8*8 + 4*4 + 4*4]

	YV12Plane prevFrame, curFrame, nextFrame;
	int fn; //frame number
	int nbx, nby;
	VecMatrix vectorsP, vectorsN;
	std::vector< std::vector<bool> > haveMVp, haveMVn; //prev, next
	std::vector< std::vector<int> > motion;
};

////////////////////////////////////////////////////////////////////////////
template<int W>
void Plane::readBlock_f(Vec v, MonoBlock<W, float> &block)
{
	const int x0 = v.x, y0 = v.y;
	assert(x0 >= -border); assert(y0 >= -border);
	assert(x0 + W <= width + border); assert(y0 + W <= height + border);

	int di = offset + y0 * stride + x0;
	int si = 0;
	const BYTE* __restrict mydata = data;
	for(int y=0; y<W; y++) {
		for(int x=0;x<W;x++)
			block.data[si + x] = mydata[di + x];
		di += stride;
		si += W;
	}
}


#endif