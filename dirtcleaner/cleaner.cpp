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
	fn = 0;
	nbx = (w + 7) / 8;
	nby = (h + 7) / 8;

	makeMatrix(vectorsP, nbx, nby);
	makeMatrix(vectorsN, nbx, nby);
	makeMatrix(haveMVp, nbx, nby);
	makeMatrix(haveMVn, nbx, nby);

	inited = 1;
}


void Cleaner::process(const VDXPixmap &src, const VDXPixmap &dst)
{
	if (fn==0) {
		prevFrame.copyFrom(src);

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

		fn++;
		return;
	}

	curFrame.copyFrom(src);
	for(int by=0;by<nby;by++)
		for(int bx=0;bx<nbx;bx++) {
			haveMVp[by][bx] = false;
			haveMVn[by][bx] = false;
		}

	auto ddata = (BYTE*)dst.data;
	const int dpitch = dst.pitch;
	auto ddata2 = (BYTE*)dst.data2;
	const int dpitch2 = dst.pitch2;
	auto ddata3 = (BYTE*)dst.data3;
	const int dpitch3 = dst.pitch3;

	const Vec d4(4,4);
	int ddd;
	//float ws[4] = {0.5,}

	for(int by=0;by<nby-1;by++) { // todo: proper border handling in last row and last column
		for(int bx=0;bx<nbx;bx++) {
			//MonoBlock<8, float> srcBlock;
			const Vec bpos(bx*8, by*8);
			//curFrame.Y.readBlock_f(v0, srcBlock);
			//Vec mv = getMVp(bx, by);
			//Vec v = bpos + mv;
			//Vec bcenter = bpos + d4;

			Vec vecs[8][8];
			for(int hy=0;hy<2;hy++) {
				for(int hx=0;hx<2;hx++) {
					Vec v0 = getMVpCenter(bx-1 + hx, by-1 + hy); // centers of blocks
					Vec v1 = getMVpCenter(bx   + hx, by-1 + hy);
					Vec v2 = getMVpCenter(bx-1 + hx, by   + hy);
					Vec v3 = getMVpCenter(bx   + hx, by   + hy);

					float kx0 = hx==0 ? 0.5 : 0.0;
					float ky0 = hy==0 ? 0.5 : 0.0;
					for(int y=0;y<4;y++) {
						float ky = y / 8.0 + ky0;
						for(int x=0;x<4;x++) {
							float kx = x / 8.0 + kx0;
							float k0 = ((1-kx)*(1-ky));
							float k1 = (kx*(1-ky));
							float k2 = ((1-kx)*ky);
							float k3 = kx*ky;
							float kk = k0+k1+k2+k3;
							assert(kk==1.0);

							FVec point = v0*k0 + v1*k1 +
								         v2*k2 + v3*k3;
							Vec ipoint(point.x + 0.5, point.y + 0.5);
							vecs[hy*4+y][hx*4+x] = ipoint;

							BYTE* p = prevFrame.Y.pixelPtr(ipoint.y, ipoint.x);
							ddata[(by*8+hy*4+y) * dpitch + bx*8 +hx*4 + x] = *p;
						}
					}
				}
			}

			/*for(int y=0;y<8;y++) {
				BYTE *p = prevFrame.Y.pixelPtr(v.y + y, v.x);
				for(int x=0;x<8;x++) {
					ddata[(by*8+y) * dpitch + bx*8 + x] = p[x];
				}
			}*/

			//Vec v = vecs[0][0];
			//Vec uv(v.x/2, v.y/2);
			//if (bx==20 && by==20)
				//ddd = 4;
				
			for(int y=0;y<4;y++) {				
				for(int x=0;x<4;x++) {
					Vec uv = vecs[y*2][x*2].div2();
					BYTE *p = prevFrame.U.pixelPtr(uv.y, uv.x);
					ddata2[(by*4+y) * dpitch2 + bx*4 + x] = *p;
					BYTE *p3 = prevFrame.V.pixelPtr(uv.y, uv.x);
					ddata3[(by*4+y) * dpitch3 + bx*4 + x] = *p3;
				}
			}//y
		}//for bx
	}//for by

	prevFrame.swap(curFrame);
	fn++;
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

Vec Cleaner::getMVpCenter(int bx, int by)
{
	return getMVp(bx, by) + Vec(bx*8 + 4, by*8 + 4);
}