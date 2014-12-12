#include "stdafx.h"
#include "cleaner.h"
#include <stdio.h>
#include <malloc.h>
#include <dvec.h>

void log(const char* str) {
	FILE *f = fopen("c:\\temp\\dc.log", "at");
	fprintf(f, "%s\n", str);
	fclose(f);
}

void logn(const char* str, int n) {
	FILE *f = fopen("c:\\temp\\dc.log", "at");
	fprintf(f, "%s = %d\n", str, n);
	fclose(f);
}

void Plane::destroy()
{
	if (ownData && data) {
		_aligned_free(data);
		data = NULL;
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
}

void Plane::copyFrom(BYTE* srcdata, int w, int h, int pitch)
{
	assert(w==width); assert(h==height);
	int si = 0;
	for(int y=0;y<h;y++) {
		memcpy( &data[offset + y * stride], &srcdata[si], w);
		si += pitch;
	}
	makeBorder(0, height);
}

void Plane::swap(Plane &other)
{
	assert(width == other.width);
	assert(height == other.height);
	assert(offset == other.offset);
	BYTE *tmp = other.data;
	other.data = data;
	data = tmp;
	bool tmpOwn = other.ownData;
	other.ownData = ownData;
	ownData = tmpOwn;
}

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
}

void YV12Plane::copyFrom(const VDXPixmap &src)
{
	Y.copyFrom((BYTE*)src.data, src.w, src.h, src.pitch);
	U.copyFrom((BYTE*)src.data2, src.w/2, src.h/2, src.pitch2);
	V.copyFrom((BYTE*)src.data3, src.w/2, src.h/2, src.pitch3);
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
	pSquad = new CSquad(si.dwNumberOfProcessors);

	inited = 1;
}

Cleaner::~Cleaner()
{
	if (pSquad) {
		delete pSquad;
		pSquad = NULL;
	}
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
	Vec vecs[8][8];
	YV12Plane &frame = prev ? prevFrame : nextFrame;
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
				for(int x=0;x<4;x++) {
					const float kx = x / 8.0 + kx0;
					const float k0 = ((1-kx)*(1-ky));
					const float k1 = (kx*(1-ky));
					const float k2 = ((1-kx)*ky);
					const float k3 = kx*ky;
					//const float kk = k0+k1+k2+k3;
					//assert(kk==1.0);

					FVec point = v0*k0 + v1*k1 +
									v2*k2 + v3*k3;
					Vec ipoint(point.x + 0.5, point.y + 0.5);
					vecs[hy*4+y][hx*4+x] = ipoint;

					BYTE* p = frame.Y.pixelPtr(ipoint.y, ipoint.x);
					yv12block[(hy*4+y) * 8 + hx*4 + x] = *p;
				}//x
			}//y
		}//hx
	}//hy

	for(int y=0;y<4;y++) {				
		for(int x=0;x<4;x++) {
			Vec uv = vecs[y*2][x*2].div2();
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

void Cleaner::RunCommand(int command, void *params, CSquadWorker *sqworker)
{
	ProcessParams *ps = (ProcessParams*) params;
	if (command==1)
		processPart(ps->src, ps->dst, ps->nFrame, sqworker);
	else 
		degrainYV12Plane(curFrame, *ps->dst, sqworker);	
}

void Cleaner::process(const VDXPixmap *src, const VDXPixmap *dst, int nFrame)
{
	ProcessParams params;
	params.src = src;
	params.dst = dst;
	params.nFrame = nFrame;
	degrainInstead = false;

	if (nFrame < 2) {		
		if (nFrame==0) {
			prevFrame.copyFrom(*src);
			copyFrame(*src, *dst);
			return;
		} else {
			curFrame.copyFrom(*src);
			degrainInstead = true;
		}		
	} else
		pSquad->RunParallel(1, &params, this);

	if (degrainInstead) 
		pSquad->RunParallel(2, &params, this);

	if (nFrame >= 2) {
		curFrame.swap(prevFrame);
		nextFrame.swap(curFrame);
	}	
}


#define PTHRESHOLD 6
#define NOISE_AMPL 25
#define NOISY 4
#define GMTHRESHOLD 0.4

void Cleaner::processPart(const VDXPixmap *pSrc, const VDXPixmap *pDst, int nFrame, CSquadWorker *sqworker)
{
	const VDXPixmap &src = *pSrc;
	const VDXPixmap &dst = *pDst;
	if (sqworker->MyNum() != 0) return;
	
	const int X = curFrame.Y.width;
	const int Y = curFrame.Y.height;
	const bool haveRightEdge = X > nbx*8;
	const bool haveLowerEdge = Y > nby*8;
	// here nFrame >= 2
	
	if (sqworker->MyNum() == 0)
		nextFrame.copyFrom(src);

	if (sqworker->MyNum() == sqworker->NumThreads()-1) {
		for(int by=0;by<nby;by++)
			for(int bx=0;bx<nbx;bx++) {
				haveMVp[by][bx] = false;
				haveMVn[by][bx] = false;
				motion[by][bx] = 0;
			}
	}

	auto ddata = (BYTE*)dst.data;
	const int dpitch = dst.pitch;
	auto ddata2 = (BYTE*)dst.data2;
	const int dpitch2 = dst.pitch2;
	auto ddata3 = (BYTE*)dst.data3;
	const int dpitch3 = dst.pitch3;

	int by0 = 0, bys = nby;
	sqworker->GetSegment(nby, by0, bys);
	const int by1 = by0 + bys;

	for(int by=by0;by<by1;by++) { 
		for(int bx=0;bx<nbx;bx++) {
			__declspec(align(16)) BYTE yv12blockP[96];
			__declspec(align(16)) BYTE yv12blockN[96];
			flowBlock(bx, by, true,  yv12blockP); 
			flowBlock(bx, by, false, yv12blockN); 

			//int uc = 128, vc = 128; //black
			int noisypixels = 0;

			for(int i=0;i<64;i++) {
				int diff = abs(yv12blockP[i] - yv12blockN[i]);
				if (diff > NOISE_AMPL) noisypixels++;
			}
			bool different = noisypixels >= NOISY;
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

		}//for bx

		if (haveRightEdge) {
			for(int y=0;y<8;y++) {
				int di = (by*8+y)*dpitch;
				BYTE *psrc = curFrame.Y.pixelPtr(by*8 + y, 0);
				for(int x=nbx*8; x<X; x++)
					ddata[di+x] = psrc[x];
			}
		}

		//simple edge check, no far propagation
		for(int bx=0; bx<nbx;bx++) {
			bool copyOrg = false;
			if (motion[by][bx]==0) {
				//left
				if (bx > 0 && motion[by][bx-1] == 1) {
					int orgDiff = vertEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = vertEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;						
						goto copyTheBlock;
					}
				}
				//right
				if (bx < nbx-1 && motion[by][bx+1] == 1) {
					int orgDiff = vertEdge(curFrame.Y.pixelPtr(by*8, (bx+1)*8), curFrame.Y.stride);
					int fltDiff = vertEdge(&ddata[by*8*dpitch + (bx+1)*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;						
						goto copyTheBlock;
					}
				}
				//up
				if (by > by0 && motion[by-1][bx] == 1) {
					int orgDiff = horEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = horEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						copyOrg = true;												
					}
				}

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
				//up
				if (by > by0 && motion[by-1][bx]==0) {
					int orgDiff = horEdge(curFrame.Y.pixelPtr(by*8, bx*8), curFrame.Y.stride);
					int fltDiff = horEdge(&ddata[by*8*dpitch + bx*8], dpitch);
					if (fltDiff > orgDiff + PTHRESHOLD) {
						motion[by-1][bx] = 1;
						//changed[by-1][bx] = 1;
						for(int y=0;y<8;y++) {
							BYTE *psrc = curFrame.Y.pixelPtr((by-1)*8 + y, bx*8);
							for(int x=0;x<8;x++) {
								ddata[((by-1)*8+y)*dpitch + bx*8 + x] = psrc[x];
							}
						}
						for(int y=0;y<4;y++) {
							BYTE *psrcU = curFrame.U.pixelPtr((by-1)*4 + y, bx*4);
							BYTE *psrcV = curFrame.V.pixelPtr((by-1)*4 + y, bx*4);
							for(int x=0;x<4;x++) {
								ddata2[((by-1)*4+y)*dpitch2 + bx*4 + x] = psrcU[x];
								ddata3[((by-1)*4+y)*dpitch3 + bx*4 + x] = psrcV[x];
							}//x
						}//y
					} // if < PTHRE
				}
			}
		}//for bx

		//copy UV for row blocks where motion==1
		for(int bx=0;bx<nbx;bx++) 
			if (motion[by][bx]==1) {
				for(int y=0;y<4;y++) {
					BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, bx*4);
					BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, bx*4);
					for(int x=0;x<4;x++) {
						ddata2[(by*4+y)*dpitch2 + bx*4 + x] = psrcU[x];
						ddata3[(by*4+y)*dpitch3 + bx*4 + x] = psrcV[x];
					}//x
				}//y
			}

		if (haveRightEdge) {
			for(int y=0;y<4;y++) {
				BYTE *psrcU = curFrame.U.pixelPtr(by*4 + y, 0);
				BYTE *psrcV = curFrame.V.pixelPtr(by*4 + y, 0);
				for(int x=nbx*4;x<X/2;x++) {
					ddata2[(by*4+y)*dpitch2 + x] = psrcU[x];
					ddata3[(by*4+y)*dpitch3 + x] = psrcV[x];
				}//x
			}//y
		}

	}//for by

	sqworker->Sync();
	if (sqworker->MyNum() != sqworker->NumThreads() - 1) return; // the rest does one thread

	int totalChanged = 0;
	for(int by=0;by<nby;by++)
		for(int bx=0;bx<nbx;bx++) {
			if (motion[by][bx]==1) { 
				totalChanged++;
			}
		}
	
	if (totalChanged >= nbx * nby * GMTHRESHOLD) {
		degrainInstead = true;
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

Vec lookAroundLP_f(Plane &plane, MonoBlock<8,float> &sblock, Vec v0, float &bestd)
{
	MonoBlock<8,float> block;
	int bestn = -1;
	for(int n=0; n<8; n++) {	
		plane.readBlock_f(v0 + around[n], block);
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


Vec findBlockHex_f(Plane &plane, MonoBlock<8,float> &sblock, Vec v0, float &bestd, bool skip0)
{
	MonoBlock<8,float> block;
	Vec cv = v0;
	if (!skip0) { // if skip0 bestd already has diff for v0
		plane.readBlock_f(v0, block);
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
			plane.readBlock_f(cv + hecks[rn], block);
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

void testStartVec(Plane &plane, MonoBlock<8,float> &srcLuma, Vec pos, Vec candidate, Vec &startVec, float &bestd)
{
	MonoBlock<8, float> block;
	Vec v = pos + candidate;
	const int W = 8, X = plane.width, Y = plane.height;
	if (v.x < -W || v.y < -W || v.x >= X || v.y >= Y) return; //out of image, not interested
	plane.readBlock_f(v, block);
	float d = fdiff(block, srcLuma);
	if (d < bestd) {
		bestd = d;
		startVec = v;
	}
}

void testCandidates(Plane &plane, MonoBlock<8,float> &srcLuma, Vec pos, int bx, int by, VecMatrix &vectors, Vec &startVec, float &bestd)
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


Vec motionSearch(int bx, int by, MonoBlock<8, float> &srcBlock, Plane &plane, VecMatrix &vectors) // => offset vector
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
	curFrame.Y.readBlock_f(v0, srcBlock);
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
	curFrame.Y.readBlock_f(v0, srcBlock);
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