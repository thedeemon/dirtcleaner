#ifndef _CLEANER_H_
#define _CLEANER_H_

#include "stdafx.h"
#include "squad.h"
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
	MonoBlock() { /*memset(data, 255, W*W);*/ } //for debug
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
	template<int W> void readBlock_f(Vec v, MonoBlock<W, float> &block);
	void swap(Plane &other);
	void makeBorder(int startRow, int nRows);
	BYTE* pixelPtr(int y, int x);

	BYTE *data; //16 bytes aligned, includes border
	int height;
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

struct ProcessParams {
	const VDXPixmap *src, *dst;
	int nFrame;
};

class Cleaner : public ISquadJob {
public:
	Cleaner() : nbx(0), nby(0), inited(0), pSquad(NULL) {}
	virtual ~Cleaner();

	void init(int w, int h);
	void process(const VDXPixmap *src, const VDXPixmap *dst, int nFrame);
	void processPart(const VDXPixmap *pSrc, const VDXPixmap *pDst, int nFrame, CSquadWorker *sqworker);

	virtual void RunCommand(int command, void *params, CSquadWorker *sqworker);


	int inited;

protected:
	Vec getMVp(int bx, int by);
	Vec getMVn(int bx, int by);
	Vec getMVCenter(int bx, int by, bool prev);
	void flowBlock(int bx, int by, bool prev, BYTE* yv12block); // yv12block [8*8 + 4*4 + 4*4]

	YV12Plane prevFrame, curFrame, nextFrame;
	int nbx, nby;
	VecMatrix vectorsP, vectorsN;
	std::vector< std::vector<bool> > haveMVp, haveMVn; //prev, next
	std::vector< std::vector<int> > motion;
	CSquad *pSquad;
	bool degrainInstead;
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