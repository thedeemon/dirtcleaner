#ifndef _CLEANER_H_
#define _CLEANER_H_

#include "stdafx.h"
#include "squad.h"
#include <stdio.h>
#include <assert.h>
#include <vd2/plugin/vdvideofilt.h>
#include <vector>

//options
#define OLD_VEC
//#define DOLOG
//#define USE_FPLANE

//algorithm parameters
#define PTHRESHOLD 6
#define NOISE_AMPL 25
#define NOISY 4
#define GMTHRESHOLD 0.4

//worker threads commands
#define CMD_PROCESS_PART 1
#define CMD_DEGRAIN      2
#define CMD_COPY         3

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

template<class T, bool Zero = true>
class TVec
{
public:
	T x, y;

	TVec(T vx, T vy) : x(vx), y(vy) {}
	TVec() { if (Zero) { x=0; y=0; } }	
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

class Plane {
public: 
	Plane() : stride(0), width(0), height(0), border(32), data(NULL), offset(0), ownData(false) {}
	virtual ~Plane();
	void create(int w, int h);
	void destroy();
	void copyFrom(BYTE* srcdata, int w, int h, int pitch, int y0, int ys);
#ifndef USE_FPLANE
	void readBlockf(Vec v, MonoBlock<8, float> &block);
#endif
	void swap(Plane &other);
	void makeBorder(int startRow, int nRows);
	BYTE* pixelPtr(int y, int x);

	BYTE *data; //16 bytes aligned, includes border
	int stride, border, offset;
	//BYTE **rows;
	int height;
	int width;
	bool ownData;
};

#ifdef USE_FPLANE
class FPlane : public Plane {
	float *fdata;
public:
	FPlane() : Plane(), fdata(NULL) {}
	~FPlane();
	void create(int w, int h);	
	void swap(FPlane &other);
	void copyFrom(BYTE* srcdata, int w, int h, int pitch, int y0, int ys);
	void readBlockf(Vec v, MonoBlock<8, float> &block);
};
#else
#define FPlane Plane
#endif

class YV12Plane {
public:
	FPlane Y;
	Plane U,V;
	int nFrame;
	void create(int w, int h);
	void destroy();
	YV12Plane() : nFrame(-11) {}
	~YV12Plane();
	void swap(YV12Plane &other);
	void copyFrom(const VDXPixmap &src, int frameNumber, int y0, int ys);
};

typedef std::vector< std::vector< Vec > > VecMatrix;

struct ProcessParams {
	const VDXPixmap *dst;
};

class Cleaner : public ISquadJob {
public:
	Cleaner() : nbx(0), nby(0), inited(0), pSquad(NULL), nThreads(0), 
		noiseAmplitude(NOISE_AMPL), maxNoisyPixels(NOISY), markUnfiltered(false) {}
	virtual ~Cleaner();

	void init(int w, int h);
	void process(VDXFBitmap *const * srcFrames, const VDXPixmap *dst, int nFrame);
	void processPart(const VDXPixmap *pDst, CSquadWorker *sqworker);

	virtual void RunCommand(int command, void *params, CSquadWorker *sqworker);

	int noiseAmplitude, maxNoisyPixels;
	bool markUnfiltered;

protected:
	Vec getMVp(int bx, int by);
	Vec getMVn(int bx, int by);
	Vec getMVCenter(int bx, int by, bool prev);
	void flowBlock(int bx, int by, bool prev, BYTE* yv12block); // yv12block [8*8 + 4*4 + 4*4]
	void parallelCopy(YV12Plane *yv12plane, const VDXPixmap* pixmap, int frameNumber);

	int nbx, nby;
	YV12Plane prevFrame, curFrame, nextFrame;
	std::vector<HANDLE> sem; //semaphores
	VecMatrix vectorsP, vectorsN;
	std::vector< std::vector<bool> > haveMVp, haveMVn; //prev, next
	std::vector< std::vector<int> > motion;
	CSquad *pSquad;
	int nThreads;
	bool degrainInstead;
public:
	int inited;
};

#endif //H