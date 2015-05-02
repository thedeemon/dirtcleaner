#include "stdafx.h"
#include "cleaner.h"
#include <stdio.h>
#include <malloc.h>
#include <dvec.h>

#include <string>

void log(const char* str) {
#ifdef DOLOG
	FILE *f = fopen("c:\\temp\\dc.log", "at");
	fprintf(f, "%s\n", str);
	fclose(f);
#endif
}

void logn(const char* str, int n) {
#ifdef DOLOG
	FILE *f = fopen("c:\\temp\\dc.log", "at");
	fprintf(f, "%s = %d\n", str, n);
	fclose(f);
#endif
}

void Plane::destroy()
{
	if (ownData && data) {
		_aligned_free(data);
		data = NULL;
		//free(rows);
		//rows = NULL;
	}
}

Plane::~Plane()
{
	destroy();
}

void Plane::create(int w, int h)
{
	destroy();
	width = w; height = h;
	stride = (w + 2*border + 15) & (~15);
	data = (BYTE*)_aligned_malloc(stride * (h + 2*border), 16);
	offset = border * stride + border;
	//rows = (BYTE**)calloc(h + 2*border, sizeof(BYTE*));
	//for(int y=-border; y < h + border; y++)
	//	rows[y + border] = &data[offset + y * stride];
	ownData = true;
}

void Plane::makeBorder(int startRow, int nRows)
{
	BYTE* const mydata = data;
	for(int y=startRow; y<startRow + nRows; y++) {
		int di = offset + y * stride;
		for(int x=0;x<border;x++)
			mydata[di - 1 - x] = mydata[di + x];
		const int n = stride - (width + border);
		for(int x=0;x<n;x++)
			mydata[di + width + x] = mydata[di + width - 1 - x];
	}
	if (startRow==0)
		for(int y=0;y<border;y++)
			memcpy(&mydata[(border-1-y)*stride], &mydata[(border+y)*stride], stride);
	if (startRow + nRows == height)
		for(int y=0;y<border;y++)
			memcpy(&mydata[(border+y+height)*stride], &mydata[(border + height - 1 -y)*stride], stride);
}

BYTE* Plane::pixelPtr(int y, int x)
{
	assert(x >= -border); assert(y >= -border);
	assert(x + 16 <= width + border); assert(y < height + border);
	return &data[offset + y * stride + x];
	//return rows[y + border] + x;
}

void Plane::copyFrom(BYTE* srcdata, int w, int h, int pitch, int y0, int ys)
{
	assert(w==width); assert(h==height);
	for(int y=y0; y<y0+ys; y++) {
		memcpy( &data[offset + y * stride], &srcdata[y * pitch], w);		
	}
	makeBorder(y0, ys);
}

void Plane::swap(Plane &other)
{
	assert(width == other.width);
	assert(height == other.height);
	assert(offset == other.offset);
	BYTE *tmp = other.data;
	other.data = data;
	data = tmp;

	//BYTE** rs = other.rows;
	//other.rows = rows;
	//rows = rs;

	bool tmpOwn = other.ownData;
	other.ownData = ownData;
	ownData = tmpOwn;
}

#ifndef USE_FPLANE
void Plane::readBlockf(Vec v, MonoBlock<8, float> &block)
{
	const int W = 8;
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

//////////////////////////////////////////////////////////////

#ifdef USE_FPLANE
FPlane::~FPlane() {
	if (fdata) {
		_aligned_free(fdata);
		fdata = NULL;
	}
}

void FPlane::create(int w, int h) {
	Plane::create(w, h);
	fdata = (float*)_aligned_malloc(stride * (h + 2*border) * sizeof(float), 16);
}

void FPlane::swap(FPlane &other) {
	Plane::swap(other);
	float *tmp = other.fdata;
	other.fdata = fdata;
	fdata = tmp;
}

void FPlane::copyFrom(BYTE* srcdata, int w, int h, int pitch, int y0, int ys) {
	assert(w==width); assert(h==height);
	for(int y=y0; y<y0+ys; y++) {
		int si = pitch * y;
		memcpy( &data[offset + y * stride], &srcdata[si], w);
		makeBorder(y, 1);
		const int di = (y + border) * stride;
		for(int x=0; x < width + border * 2; x++)
			fdata[di + x] = data[di + x];
	}
	if (y0==0) {
		for(int x=0; x < border * stride; x++)
			fdata[x] = data[x];
	}
	if (y0 + ys == height) {
		const int di = (border + height) * stride;
		for(int x=0; x < border * stride; x++)
			fdata[di + x] = data[di + x];
	}
}

void FPlane::readBlockf(Vec v, MonoBlock<8, float> &block)
{
	const int W = 8;
	const int x0 = v.x, y0 = v.y;
	assert(x0 >= -border); assert(y0 >= -border);
	assert(x0 + W <= width + border); assert(y0 + W <= height + border);

	int di = offset + y0 * stride + x0;
	int si = 0;
	const float* __restrict mydata = fdata;
	for(int y=0; y<W; y++) {
		for(int x=0;x<W;x++)
			block.data[si + x] = mydata[di + x];
		di += stride;
		si += W;
	}
}
#endif

//////////////////////////////////////////////////////////////

void YV12Plane::create(int w, int h) 
{
	Y.create(w,h);
	U.create(w/2, h/2);
	V.create(w/2, h/2);
}

void YV12Plane::destroy()
{
	Y.destroy(); U.destroy(); V.destroy();
}

YV12Plane::~YV12Plane()
{
	destroy();
}

void YV12Plane::swap(YV12Plane &other)
{
	Y.swap(other.Y);	U.swap(other.U);	V.swap(other.V);
	int t = other.nFrame;
	other.nFrame = nFrame;
	nFrame = t;
}

void YV12Plane::copyFrom(const VDXPixmap &src, int frameNumber, int y0, int ys)
{
	assert((y0 & 1)==0); assert((ys & 1)==0); assert(ys > 1);
	Y.copyFrom((BYTE*)src.data, src.w, src.h, src.pitch, y0, ys);
	U.copyFrom((BYTE*)src.data2, src.w/2, src.h/2, src.pitch2, y0/2, ys/2);
	V.copyFrom((BYTE*)src.data3, src.w/2, src.h/2, src.pitch3, y0/2, ys/2);
	nFrame = frameNumber;
}

template<class T>
void makeMatrix(std::vector< std::vector<T> > &mat, int nbx, int nby)
{
	mat.resize(nby);
	for(int y=0;y<nby;y++) {
		mat[y].clear();
		mat[y].resize(nbx);
	}
}

void Cleaner::init(int w, int h) 
{
	prevFrame.create(w,h); curFrame.create(w,h); nextFrame.create(w,h);
	nbx = w / 8;
	nby = h / 8;

	makeMatrix(vectorsP, nbx, nby);
	makeMatrix(vectorsN, nbx, nby);
	makeMatrix(haveMVp, nbx, nby);
	makeMatrix(haveMVn, nbx, nby);
	makeMatrix(motion, nbx, nby);

	SYSTEM_INFO si;
	GetSystemInfo(&si);
	nThreads = si.dwNumberOfProcessors;
	pSquad = new CSquad(nThreads);
	
	sem.resize(nThreads);
	for(int i=0;i<nThreads;i++)
		sem[i] = CreateSemaphore(NULL,0,nbx,NULL);

	inited = 1;
}

Cleaner::~Cleaner()
{
	if (pSquad) {
		delete pSquad;
		pSquad = NULL;
	}
	for(int i=0;i<sem.size();i++)
		CloseHandle(sem[i]);
}

void copyFrame(const VDXPixmap &src, const VDXPixmap &dst)
{
	{
		auto ddata = (BYTE*)dst.data;
		auto sdata = (BYTE*)src.data;
		for(int y=0;y<dst.h;y++) {
			//for(int x=0;x<dst.w;x++)
			//	ddata[y * dst.pitch + x] = x ^ y;
			memcpy(&ddata[y * dst.pitch], &sdata[y * src.pitch], dst.w);
		}
	}
	{
		auto ddata2 = (BYTE*)dst.data2;
		auto sdata2 = (BYTE*)src.data2;
		for(int y=0;y<dst.h/2;y++) 
			memcpy(&ddata2[y * dst.pitch2], &sdata2[y * src.pitch2], dst.w/2);
	}
	{
		auto ddata3 = (BYTE*)dst.data3;
		auto sdata3 = (BYTE*)src.data3;
		for(int y=0;y<dst.h/2;y++) 
			memcpy(&ddata3[y * dst.pitch3], &sdata3[y * src.pitch3], dst.w/2);
	}
}

#define CMP(x,y) if (x > y) { int t = x; x = y;  y = t; }

void degrainPlane(BYTE *src, int spitch, BYTE *dst, int dpitch, int X, int Y)
{
	memcpy(dst, src, X);
	for(int y=1; y<Y-1; y++) {
		dst[y*dpitch] = src[y*spitch];
		for(int x=1;x < X-1;x++) {
			int a[8];
			int si = y * spitch + x;
			a[0] = src[si-spitch-1];
			a[1] = src[si-spitch];
			a[2] = src[si-spitch+1];

			a[3] = src[si-1];
			int c = src[si];
			a[4] = src[si+1];

			a[5] = src[si+spitch-1];
			a[6] = src[si+spitch];
			a[7] = src[si+spitch+1];

			CMP(a[0], a[1]); CMP(a[2], a[3]); CMP(a[4], a[5]); CMP(a[6], a[7]); //16 comparisons
			CMP(a[0], a[2]); CMP(a[1], a[3]); CMP(a[4], a[6]); CMP(a[5], a[7]);
			CMP(a[1], a[2]); CMP(a[5], a[6]); CMP(a[0], a[4]); CMP(a[3], a[7]);
			CMP(a[1], a[5]); CMP(a[2], a[6]);
			CMP(a[1], a[4]); CMP(a[3], a[6]);

			dst[y * dpitch + x] = max(min(c, a[6]), a[1]);
		}
		dst[y*dpitch + X-1] = src[y*spitch + X-1];
	}
	memcpy(&dst[(Y-1)*dpitch], &src[(Y-1)*spitch], X);
}

void degrainFrame(const VDXPixmap &src, const VDXPixmap &dst)
{
	degrainPlane((BYTE*)src.data,  src.pitch,  (BYTE*)dst.data,  dst.pitch,  dst.w,   dst.h);
	degrainPlane((BYTE*)src.data2, src.pitch2, (BYTE*)dst.data2, dst.pitch2, dst.w/2, dst.h/2);
	degrainPlane((BYTE*)src.data3, src.pitch3, (BYTE*)dst.data3, dst.pitch3, dst.w/2, dst.h/2);
}

void degrainYV12Plane(YV12Plane &src, const VDXPixmap &dst, CSquadWorker *sqworker)
{
	const int mynum = sqworker->MyNum();
	if (mynum == 0)
		degrainPlane(src.Y.pixelPtr(0,0),  src.Y.stride, (BYTE*)dst.data,  dst.pitch,  dst.w,   dst.h);
	if (mynum == sqworker->NumThreads() - 1) {
		degrainPlane(src.U.pixelPtr(0,0),  src.U.stride, (BYTE*)dst.data2, dst.pitch2, dst.w/2, dst.h/2);
		degrainPlane(src.V.pixelPtr(0,0),  src.V.stride, (BYTE*)dst.data3, dst.pitch3, dst.w/2, dst.h/2);
	}
}

void Cleaner::flowBlock(int bx, int by, bool prev, BYTE* yv12block) // yv12block [8*8 + 4*4 + 4*4]
{
	TVec<int, false> vecs[8][8];
	//Vec vecs2[8][8];
	YV12Plane &frame = prev ? prevFrame : nextFrame;
	const float dkx = 0.125;

	for(int hy=0;hy<2;hy++) {
		const float ky0 = hy==0 ? 0.5 : 0.0;
		for(int hx=0;hx<2;hx++) {
			const Vec v0 = getMVCenter(bx-1 + hx, by-1 + hy, prev); // centers of blocks
			const Vec v1 = getMVCenter(bx   + hx, by-1 + hy, prev);
			const Vec v2 = getMVCenter(bx-1 + hx, by   + hy, prev);
			const Vec v3 = getMVCenter(bx   + hx, by   + hy, prev);

			const float kx0 = hx==0 ? 0.5 : 0.0;			
			for(int y=0;y<4;y++) {
				const float ky = y / 8.0 + ky0;
#ifndef OLD_VEC
				const float dk0 = -dkx * (1 - ky);
				const float dk1 = dkx * (1 - ky);
				const float dk2 = -dkx * ky;
				const float dk3 = dkx * ky;
				const FVec dv = v0 * dk0 + v1 * dk1 + v2 * dk2 + v3 * dk3;
				const FVec point2 = v0*((1-kx0)*(1-ky)) + v1*(kx0*(1-ky)) + v2 * ((1-kx0)*ky) + v3 * (kx0*ky);
#endif
				for(int x=0;x<4;x++) {
#ifdef OLD_VEC
					const float kx = x / 8.0 + kx0;
					const float k0 = (1-kx)*(1-ky);
					const float k1 = kx*(1-ky);
					const float k2 = (1-kx)*ky;
					const float k3 = kx*ky;
					//const float kk = k0+k1+k2+k3;
					//assert(kk==1.0);

					FVec point = v0*k0 + v1*k1 + v2*k2 + v3*k3;
					

					//FVec point2 = v0*k0_ + v1*k1_ + v2*k2_ + v3*k3_ = point + dv
					// dv = v0 * (k0_ - k0) + v1 * (k1_ - k1) + v2 * (k2_ - k2) + v3 * (k3_ - k3)
					// dkx = kx_ - kx = 0.125
					// dk0 = (1-kx_)*(1-ky) - (1-kx)*(1-ky) = (1-kx_ -1 + kx) * (1-ky) = -dkx * (1 - ky)
					// dk1 = kx_*(1-ky) - kx*(1-ky) = dkx * (1 - ky)
					// dk2 = (1-kx_)*ky - (1-kx)*ky = -dkx * ky
					// dk3 = dkx * ky
#else
					FVec point = point2 + dv * x;
#endif
					TVec<int, false> ipoint(point.x + 0.5, point.y + 0.5);
					vecs[hy*4+y][hx*4+x] = ipoint;

					//assert(ipoint2 == ipoint);
					BYTE* p = frame.Y.pixelPtr(ipoint.y, ipoint.x);
					yv12block[(hy*4+y) * 8 + hx*4 + x] = *p;
					//point2 += dv;
				}//x
			}//y
		}//hx
	}//hy

	for(int y=0;y<4;y++) {				
		for(int x=0;x<4;x++) {
			TVec<int, false> uv = vecs[y*2][x*2].div2();
			BYTE *p = frame.U.pixelPtr(uv.y, uv.x);
			yv12block[64 + y * 4 + x] = *p;
			BYTE *p3 = frame.V.pixelPtr(uv.y, uv.x);
			yv12block[80 + y * 4 + x] = *p3;
		}
	}//y
}

int vertEdge(BYTE *p, int pitch) // compares with lefter column
{
	int sum = 0;
	for(int i=0;i<8;i++) {
		sum += abs(p[0] - p[-1]);
		p += pitch;
	}
	return sum;
}

int horEdge(BYTE *p, int pitch) //compares with upper line 
{
	int sum = 0;
	BYTE *up = p - pitch;
	for(int i=0;i<8;i++)
		sum += abs(p[i] - up[i]);
	return sum;
}

struct CopyParams {
	YV12Plane *yv12plane; const VDXPixmap* pixmap; int frameNumber;
};


void Cleaner::RunCommand(int command, void *params, CSquadWorker *sqworker) {
	ProcessParams *ps = (ProcessParams*) params;
	switch(command) {
		case CMD_PROCESS_PART: 	processPart(ps->dst, sqworker); break;
		case CMD_DEGRAIN: 		degrainYV12Plane(curFrame, *ps->dst, sqworker);	break;
		case CMD_COPY: {
			CopyParams *cp = (CopyParams*) params;
			int H = cp->yv12plane->Y.height;
			int y0 = 0, ys = H;
			sqworker->GetSegment(H/2, y0, ys);
			y0 *= 2; ys *= 2;
			if (ys > 1)
				cp->yv12plane->copyFrom(*cp->pixmap, cp->frameNumber, y0, ys);
		}
	}
}

void Cleaner::parallelCopy(YV12Plane *yv12plane, const VDXPixmap* pixmap, int frameNumber) {
	CopyParams ps;
	ps.yv12plane = yv12plane; ps.pixmap = pixmap; ps.frameNumber = frameNumber;
	pSquad->RunParallel(CMD_COPY, &ps, this);
}

/*void calcDiffHist(Plane &a, Plane &b, int* tab) { //tab[256]
	assert(a.height == b.height); assert(a.width == b.width);
	for(int y=0;y<a.height;y++) {
		const BYTE * const pa = a.pixelPtr(y,0);
		const BYTE * const pb = b.pixelPtr(y,0);
		for(int x=0;x<a.width;x++) {
			int d = abs(pa[x] - pb[x]); 
			tab[d]++;
		}
	}
}

void frameDiffHist(YV12Plane &a, YV12Plane &b, int *ytab, int *uvtab) {
	calcDiffHist(a.Y, b.Y, ytab);
	calcDiffHist(a.U, b.U, uvtab);
	calcDiffHist(a.V, b.V, uvtab);
}*/

bool samePlane(Plane &a, Plane &b) {
	assert(a.height == b.height); assert(a.width == b.width);
	for(int y=0;y<a.height;y++)  {
		const BYTE * const pa = a.pixelPtr(y,0);
		const BYTE * const pb = b.pixelPtr(y,0);
		for(int x=0;x<a.width;x++) {
			const int d = abs(pa[x] - pb[x]); 
			if (d > 4) return false;
		}
	}
	return true;
}

bool sameFrames(YV12Plane &a, YV12Plane &b) {
	return samePlane(a.Y, b.Y) && samePlane(a.U, b.U) && samePlane(a.V, b.V);
}

void Cleaner::process(VDXFBitmap *const *srcFrames, const VDXPixmap *dst, int nFrame) {
	ProcessParams params;
	params.dst = dst;
	degrainInstead = false;

	int df = nFrame - curFrame.nFrame;
	if (df == 1) {
		if (sameFrames(curFrame, nextFrame)) { //ABB -> ABC
			logn("cur === next for cur.nFrame", nFrame);
			curFrame.nFrame = nextFrame.nFrame;
		} else { //ABC -> BCD
			curFrame.swap(prevFrame);
			nextFrame.swap(curFrame);
			//nextFrame.copyFrom(*srcFrames[2]->mpPixmap, srcFrames[2]->mFrameNumber);					
		}
		parallelCopy(&nextFrame, srcFrames[2]->mpPixmap, srcFrames[2]->mFrameNumber);
	} else
	if (df == -1) {
		curFrame.swap(nextFrame);
		prevFrame.swap(curFrame);
		//prevFrame.copyFrom(*srcFrames[0]->mpPixmap, srcFrames[0]->mFrameNumber);	
		parallelCopy(&prevFrame, srcFrames[0]->mpPixmap, srcFrames[0]->mFrameNumber);
	} else 
	if (abs(df) > 1) { // seek happened (or just the very first frame)
		//prevFrame.copyFrom(*srcFrames[0]->mpPixmap, srcFrames[0]->mFrameNumber);		
		//curFrame.copyFrom(*srcFrames[1]->mpPixmap, srcFrames[1]->mFrameNumber);		
		//nextFrame.copyFrom(*srcFrames[2]->mpPixmap, srcFrames[2]->mFrameNumber);		
		#ifdef DOLOG
		log(" process: parallelCopy 3 input frames");
		char str[1024];
		sprintf(str, " srcFrames[0] pixmap=%X data=%X data2=%X", srcFrames[0]->mpPixmap, srcFrames[0]->mpPixmap->data, srcFrames[0]->mpPixmap->data2);
		log(str);
		sprintf(str, " srcFrames[1] pixmap=%X data=%X data2=%X", srcFrames[1]->mpPixmap, srcFrames[1]->mpPixmap->data, srcFrames[1]->mpPixmap->data2);
		log(str);
		sprintf(str, " srcFrames[2] pixmap=%X data=%X data2=%X", srcFrames[2]->mpPixmap, srcFrames[2]->mpPixmap->data, srcFrames[2]->mpPixmap->data2);
		log(str);
		#endif

		parallelCopy(&prevFrame, srcFrames[0]->mpPixmap, srcFrames[0]->mFrameNumber);
		parallelCopy(&curFrame, srcFrames[1]->mpPixmap, srcFrames[1]->mFrameNumber);
		parallelCopy(&nextFrame, srcFrames[2]->mpPixmap, srcFrames[2]->mFrameNumber);
		makeMatrix(vectorsP, nbx, nby);
		makeMatrix(vectorsN, nbx, nby);
	}

	/*int pnFrame = nFrame == 0 ? 0 : nFrame - 1;
	assert(prevFrame.nFrame == pnFrame);
	assert(curFrame.nFrame == nFrame);
	assert(nextFrame.nFrame == nFrame + 1);*/

	if (nFrame < 2) {		
		//copyFrame(*srcFrames[1]->mpPixmap, *dst);
		degrainInstead = true;
		log("nFrame < 2 => degrain");
	} else {
		for(int by=0;by<nby;by++)
			for(int bx=0;bx<nbx;bx++) {
				haveMVp[by][bx] = false;
				haveMVn[by][bx] = false;
				motion[by][bx] = 0;
			}
		
		pSquad->RunParallel(CMD_PROCESS_PART, &params, this);
	}

	if (degrainInstead) 
		pSquad->RunParallel(CMD_DEGRAIN, &params, this);
}

void Cleaner::processPart(const VDXPixmap *pDst, CSquadWorker *sqworker)
{
	//const VDXPixmap &src = *pSrc;
	const VDXPixmap &dst = *pDst;
	
	const int X = curFrame.Y.width;
	const int Y = curFrame.Y.height;
	const bool haveRightEdge = X > nbx*8;
	const bool haveLowerEdge = Y > nby*8;
	// here nFrame >= 2
	const int lastThread = sqworker->NumThreads() - 1;
	const int myNum = sqworker->MyNum();
	
	auto ddata = (BYTE*)dst.data;
	const int dpitch = dst.pitch;
	auto ddata2 = (BYTE*)dst.data2;
	const int dpitch2 = dst.pitch2;
	auto ddata3 = (BYTE*)dst.data3;
	const int dpitch3 = dst.pitch3;

	int by0 = 0, bys = nby;
	sqworker->GetSegment(nby, by0, bys);
	const int by1 = by0 + bys;

	for(int bx=0;bx<nbx;bx++) {
		if (myNum>0)
			WaitForSingleObject(sem[myNum-1], INFINITE);
		for(int by=by0;by<by1;by++) { 		
			__declspec(align(16)) BYTE yv12blockP[96];
			__declspec(align(16)) BYTE yv12blockN[96];
			flowBlock(bx, by, true,  yv12blockP); 
			flowBlock(bx, by, false, yv12blockN); 

			//int uc = 128, vc = 128; //black
			int noisypixels = 0;

			for(int i=0;i<64;i++) {
				int diff = abs(yv12blockP[i] - yv12blockN[i]);
				if (diff > noiseAmplitude) noisypixels++;
			}
			bool different = noisypixels >= maxNoisyPixels;
			motion[by][bx] = different ? 1 : 0; 

			if (!different) {//motion=0
				for(int y=0;y<8;y++) {
					BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, bx*8);
					for(int x=0;x<8;x++) {
						const int i = y*8 + x;
						int c = psrc[x];
						int upper = max(yv12blockN[i], yv12blockP[i]);
						int lower = min(yv12blockN[i], yv12blockP[i]);
						int cl = max(c, lower);
						ddata[(by*8+y)*dpitch + bx*8 + x] = min(cl, upper);
					}
				}
				for(int y=0;y<4;y++) {
					BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, bx*4);
					BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, bx*4);
					for(int x=0;x<4;x++) {
						const int j = y*4 + x;
						{
							int i = j + 64;
							int upper = max(yv12blockN[i], yv12blockP[i]);
							int lower = min(yv12blockN[i], yv12blockP[i]);
							int c = psrcU[x];
							int cl = max(c, lower);
							ddata2[(by*4+y)*dpitch2 + bx*4 + x] = min(cl, upper);
						}
						{
							int i = j + 80;
							int upper = max(yv12blockN[i], yv12blockP[i]);
							int lower = min(yv12blockN[i], yv12blockP[i]);
							int c = psrcV[x];
							int cl = max(c, lower);
							ddata3[(by*4+y)*dpitch3 + bx*4 + x] = min(cl, upper);
						}
					}//x
				}//y
			} else { ////motion=1 : blocks are too different, copy source
				for(int y=0;y<8;y++) {
					BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, bx*8);
					for(int x=0;x<8;x++) {
						ddata[(by*8+y)*dpitch + bx*8 + x] = psrc[x];
					}
				}
			}//different?
		}//for by

		if (myNum < nThreads-1)
		    ReleaseSemaphore(sem[myNum], 1, NULL);

		//simple edge check, no far propagation
		for(int by=by0;by<by1;by++) {
			bool copyOrg = false;
			if (motion[by][bx]==0) {
				//up
				if (by > by0 && motion[by-1][bx] == 1) {
					int orgDiff = horEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = horEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;												
						goto copyTheBlock;
					}
				}
				//down
				if (by < by1-1 && motion[by+1][bx] == 1) {
					int orgDiff = horEdge(curFrame.Y.pixelPtr((by+1)*8, bx*8), curFrame.Y.stride);
					int fltDiff = horEdge(&ddata[(by+1)*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;		
						goto copyTheBlock;
					}
				}
				//left
				if (bx > 0 && motion[by][bx-1] == 1) {
					int orgDiff = vertEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = vertEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;						
						goto copyTheBlock;
					}
				}
				//right block is not calculated yet
				/*if (bx < nbx-1 && motion[by][bx+1] == 1) {
					int orgDiff = vertEdge(curFrame.Y.pixelPtr(by*8, (bx+1)*8), curFrame.Y.stride);
					int fltDiff = vertEdge(&ddata[by*8*dpitch + (bx+1)*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;						
						goto copyTheBlock;
					}
				}*/

				if (copyOrg) {
					copyTheBlock:
					motion[by][bx] = 1;
					for(int y=0;y<8;y++) {
						BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, bx*8);
						for(int x=0;x<8;x++) {
							ddata[(by*8+y)*dpitch + bx*8 + x] = psrc[x];
						}
					}
				}
			} else { //this_block.motion=1
				//left
				if (bx > 0 && motion[by][bx-1]==0) {
					int orgDiff = vertEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = vertEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						motion[by][bx-1] = 1;
						for(int y=0;y<8;y++) {
							BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, (bx-1)*8);
							for(int x=0;x<8;x++) {
								ddata[(by*8+y)*dpitch + (bx-1)*8 + x] = psrc[x];
							}
						}
						for(int y=0;y<4;y++) {
							BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, (bx-1)*4);
							BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, (bx-1)*4);
							for(int x=0;x<4;x++) {
								ddata2[(by*4+y)*dpitch2 + (bx-1)*4 + x] = psrcU[x];
								ddata3[(by*4+y)*dpitch3 + (bx-1)*4 + x] = psrcV[x];
							}//x
						}//y
						if (markUnfiltered)
							ddata3[(by*4)*dpitch3 + (bx-1)*4] = 0; //dbg					
					} // if < PTHRE
				}
			}//if motion==1
		}//for by

		//copy UV for column blocks where motion==1
		for(int by=by0;by<by1;by++) 
			if (motion[by][bx]==1) {
				for(int y=0;y<4;y++) {
					BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, bx*4);
					BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, bx*4);
					for(int x=0;x<4;x++) {
						ddata2[(by*4+y)*dpitch2 + bx*4 + x] = psrcU[x];
						ddata3[(by*4+y)*dpitch3 + bx*4 + x] = psrcV[x];
					}//x
				}//y
				if (markUnfiltered)
					ddata3[(by*4)*dpitch3 + bx*4] = 0;//dbg				
			}
	}//for bx

	if (haveRightEdge) {
		for(int by=by0;by<by1;by++) {
			for(int y=0;y<8;y++) {
				const int di = (by*8+y)*dpitch;
				BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, 0);
				for(int x=nbx*8; x<X; x++)
					ddata[di+x] = psrc[x];
			}
			for(int y=0;y<4;y++) {
				BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, 0);
				BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, 0);
				for(int x=nbx*4;x<X/2;x++) {
					ddata2[(by*4+y)*dpitch2 + x] = psrcU[x];
					ddata3[(by*4+y)*dpitch3 + x] = psrcV[x];
				}//x
			}//y
		}//by
	}//if right edge


	sqworker->Sync();
	if (myNum != lastThread) return; // the rest does one last thread

	int totalChanged = 0;
	for(int by=0;by<nby;by++)
		for(int bx=0;bx<nbx;bx++) {
			if (motion[by][bx]==1) { 
				totalChanged++;
			}
		}
	
	SHOW(totalChanged);
	if (totalChanged >= nbx * nby * GMTHRESHOLD) {
		degrainInstead = true;
		log("degrain = true");
	} else {
		if (haveLowerEdge) {
			for(int y=nby*8; y<Y; y++)
				memcpy(&ddata[y*dpitch], curFrame.Y.pixelPtr(y, 0), X);
			for(int y=nby*4; y<Y/2; y++) {
				memcpy(&ddata2[y*dpitch2], curFrame.U.pixelPtr(y, 0), X/2);
				memcpy(&ddata3[y*dpitch3], curFrame.V.pixelPtr(y, 0), X/2);
			}
		}
	}
}


template <int W>
float fdiff(MonoBlock<W,float> &block, MonoBlock<W,float> &sblock)
{
	F32vec4 *a = (F32vec4 *)block.data;
	F32vec4 *b = (F32vec4 *)sblock.data;
	F32vec4 s(0.0);
	for(int i=0; i<W*W/4; i++) {
		const F32vec4 t = a[i] - b[i];
		s += t*t;
	}
	return add_horizontal(s);
}

static Vec around[8] =   { Vec(-1,-1), Vec(0,-1), Vec(1,-1), 
						   Vec(-1,0),             Vec(1,0), 
	                       Vec(-1,1),  Vec(0,1),  Vec(1,1) };

static Vec hecks[6] = { Vec(1,-2), Vec(2,0), Vec(1,2), 
	                    Vec(-1,2), Vec(-2,0), Vec(-1,-2) };

Vec lookAroundLP_f(FPlane &plane, MonoBlock<8,float> &sblock, Vec v0, float &bestd)
{
	MonoBlock<8,float> block;
	int bestn = -1;
	for(int n=0; n<8; n++) {	
		plane.readBlockf(v0 + around[n], block);
		auto d = fdiff(block, sblock);
		if (d<bestd) {
			bestd = d;
			bestn = n;
		}
	}
	if (bestn>=0) 
		v0.move( around[bestn] );
	return v0;
}


Vec findBlockHex_f(FPlane &plane, MonoBlock<8,float> &sblock, Vec v0, float &bestd, bool skip0)
{
	MonoBlock<8,float> block;
	Vec cv = v0;
	if (!skip0) { // if skip0 bestd already has diff for v0
		plane.readBlockf(v0, block);
		bestd = fdiff(block, sblock);
	}
	int start = 0, end = 6, lastbestn, bestn = -1; 
	const int W = 8, X = plane.width, Y = plane.height;

	while(true) {
		lastbestn = bestn;
		bestn = -1;
		if (cv.x < -W || cv.y < -W || cv.x > X || cv.y > Y) break;
		for(int n=start;n<end;n++) {
			int rn = n % 6;
			plane.readBlockf(cv + hecks[rn], block);
			auto d = fdiff(block, sblock);
			if (d < bestd) {
				bestd = d;
				bestn = rn;
			}
		}
		if (bestn >= 0) { //move
			start = bestn + 6 - 1;
			end = bestn + 6 + 2;
			cv += hecks[bestn];
		} else
			break;
	}
	cv = lookAroundLP_f(plane, sblock, cv, bestd);
	return cv;
}

void testStartVec(FPlane &plane, MonoBlock<8,float> &srcLuma, Vec pos, Vec candidate, Vec &startVec, float &bestd)
{
	MonoBlock<8, float> block;
	Vec v = pos + candidate;
	const int W = 8, X = plane.width, Y = plane.height;
	if (v.x < -W || v.y < -W || v.x >= X || v.y >= Y) return; //out of image, not interested
	plane.readBlockf(v, block);
	float d = fdiff(block, srcLuma);
	if (d < bestd) {
		bestd = d;
		startVec = v;
	}
}

void testCandidates(FPlane &plane, MonoBlock<8,float> &srcLuma, Vec pos, int bx, int by, VecMatrix &vectors, Vec &startVec, float &bestd)
{
	Vec vs[6], cands[6];
	int nc = 1;
	cands[0] = Vec(0,0);
	cands[nc++] = vectors[by][bx]; // old one
	if (by > 0) {
		cands[nc++] = vectors[by-1][bx];
		//cands[nc++] = vectors[by-1][bx+1];
	}
	if (bx > 0) {
		cands[nc++] = vectors[by][bx-1];
		//cands[nc++] = vectors[hby+1][hbx-1];
	}
	if (bx > 0 && by > 0)
		cands[nc++] = vectors[by-1][bx-1];
	vs[0] = cands[0];
	int nv = 1;
	for(int i=1;i<nc;i++) {
		bool dup = false;
		for(int j=0;j<i;j++)
			if (cands[i] == vs[j]) { dup = true; break;	}
		if (!dup)
			vs[nv++] = cands[i];
	}

	for(int i=0;i<nv;i++)
		testStartVec(plane, srcLuma, pos, vs[i], startVec, bestd);
}


Vec motionSearch(int bx, int by, MonoBlock<8, float> &srcBlock, FPlane &plane, VecMatrix &vectors) // => offset vector
{
	float bestd = 10000000;	
	const Vec v0(bx*8, by*8);
	Vec startVec = v0;	
	testCandidates(plane, srcBlock, v0, bx, by, vectors, startVec, bestd);
	Vec v1 = findBlockHex_f(plane, srcBlock, startVec, bestd, true);
	//comp_1HsQ(srcLuma, v0, v1, hbx, hby, vectors, diffs); // 1 = comp_1HQsQ seems best for smooth among 32; 8 = comp_1HsQ
	//comp(srcLuma, fsrcLuma, v0, v1, hbx, hby, vectors, diffs);*/
	Vec dv = v1 - v0;
	return dv;
}

Vec Cleaner::getMVp(int bx, int by)
{
	if (bx < 0) bx = 0;
	if (bx >= nbx) bx = nbx - 1;
	if (by < 0) by = 0;
	if (by >= nby) by = nby - 1;

	if (haveMVp[by][bx]) return vectorsP[by][bx];
	
	MonoBlock<8, float> srcBlock;
	const Vec v0(bx*8, by*8);
	curFrame.Y.readBlockf(v0, srcBlock);
	Vec mv = motionSearch(bx, by, srcBlock, prevFrame.Y, vectorsP);
	vectorsP[by][bx] = mv;
	haveMVp[by][bx] = true;
	return mv;
}

Vec Cleaner::getMVn(int bx, int by)
{
	if (bx < 0) bx = 0;
	if (bx >= nbx) bx = nbx - 1;
	if (by < 0) by = 0;
	if (by >= nby) by = nby - 1;

	if (haveMVn[by][bx]) return vectorsN[by][bx];
	
	MonoBlock<8, float> srcBlock;
	const Vec v0(bx*8, by*8);
	curFrame.Y.readBlockf(v0, srcBlock);
	Vec mv = motionSearch(bx, by, srcBlock, nextFrame.Y, vectorsN);
	vectorsN[by][bx] = mv;
	haveMVn[by][bx] = true;
	return mv;
}

Vec Cleaner::getMVCenter(int bx, int by, bool prev)
{
	Vec v = prev ? getMVp(bx, by) : getMVn(bx, by);
	return v + Vec(bx*8 + 4, by*8 + 4);
}