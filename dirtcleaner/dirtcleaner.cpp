#include "stdafx.h"
#include <vd2/plugin/vdvideofilt.h>
#include <WindowsX.h>
#include <stdio.h>
#include "cleaner.h"

struct MyFilterData {
	Cleaner *cln() {
		if (pCleaner==NULL) 
			pCleaner = new Cleaner();
		return pCleaner;
	}

	void deinit() {
		if (pCleaner) {
			delete pCleaner;
			pCleaner = NULL;
		}
	}
private:
	Cleaner *pCleaner;
};

long paramProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	const VDXPixmapLayout& pxsrc = *fa->src.mpPixmapLayout;
          VDXPixmapLayout& pxdst = *fa->dst.mpPixmapLayout;
          
    // check for a source format that we support
    if (pxsrc.format != nsVDXPixmap::kPixFormat_YUV420_Planar)
        return FILTERPARAM_NOT_SUPPORTED;

    // set old depth value to zero to indicate new pixmap layout should be used
    fa->dst.depth = 0;

	return FILTERPARAM_SWAP_BUFFERS | FILTERPARAM_SUPPORTS_ALTFORMATS | FILTERPARAM_HAS_LAG(1);
}

int runProc(const VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("runProc");
	auto src = fa->src.mpPixmap;
	auto dst = fa->dst.mpPixmap;
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	Cleaner *c = pData->cln();
	if (c->inited==0) 
		c->init(src->w, src->h);
	VDXFilterStateInfo *psi = fa->pfsi;
	int fn = 0;
	if (psi) {
		SHOW(psi->lCurrentFrame);
		fn = psi->lCurrentFrame;
	} else log("pfsi is null");
	c->process(*src, *dst, fn);	
	return 0;
}

int startProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("startProc");
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	auto src = fa->src.mpPixmap;
	Cleaner *c = pData->cln();
	c->init(src->w, src->h);
	return 0;
}

int endProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("endProc");
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	pData->deinit();
	return 0;
}

int configProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff, VDXHWND hwndParent) {
    MyFilterData* pData = (MyFilterData*)fa->filter_data;
	//pData->ifp = fa->ifp;
    //auto res = !DialogBoxParam(g_hInst, MAKEINTRESOURCE(IDD_SETTINGS), (HWND)hwndParent, SettingsDlgProc, (LPARAM)pData);
    return 0;// res;
}

static void stringProc2(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, char *buf, int maxlen) {
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	//_snprintf(buf, maxlen, " (%d-%d)", pData->targetMin, pData->targetMax );
}

void stringProc(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, char *buf) {
	stringProc2(fa, ff, buf, 80);
}

static struct VDXFilterDefinition myfilter_definition={
    0,0,NULL,                       // next, prev, and module (set to zero)
    "Film Dirt Cleaner",                    // name
    "Motion-compensated temporal dirt cleaner.",    // description
    "Infognition Co. Ltd.",         // author / maker
    NULL,                           // no private data
    sizeof(MyFilterData),           // no instance data size
    NULL,                           // no initProc
    NULL,                           // no deinitProc
    runProc,                   // runProc
	paramProc,
	configProc,
	stringProc,
	startProc,
	endProc,
	0,0, //&script_obj, fssProc,
	stringProc2,0,0,0,
	0, 0,
	0,0    
};


VDXFilterDefinition *g_registeredFilterDef;

extern "C" __declspec(dllexport) int __cdecl VirtualdubFilterModuleInit2(struct VDXFilterModule *fm, const VDXFilterFunctions *ff, int& vdfd_ver, int& vdfd_compat) 
{
	/*if (vdfd_ver < 14) {
		MessageBox(NULL, "VirtualDub version 1.9.1 or higher required.", "Infognition Super Resolution", MB_ICONERROR);
		return 1;
	}*/
    g_registeredFilterDef = ff->addFilter(fm, &myfilter_definition, sizeof(VDXFilterDefinition));

    vdfd_ver        = VIRTUALDUB_FILTERDEF_VERSION; //14
    vdfd_compat     = 14;    // VD 1.9.1 or higher

    return 0;
}

extern "C" __declspec(dllexport) void __cdecl VirtualdubFilterModuleDeinit(struct VDXFilterModule *fm, const VDXFilterFunctions *ff) 
{
    ff->removeFilter(g_registeredFilterDef);
}
