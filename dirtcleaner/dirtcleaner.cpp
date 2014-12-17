#include "stdafx.h"
#include <vd2/plugin/vdvideofilt.h>
#include <WindowsX.h>
#include <stdio.h>
#include "cleaner.h"
#include "resource.h"
#include <commctrl.h>
#include <shellapi.h>

struct MyFilterData {
	int noiseLevel, numPixels;
	bool markUnfiltered;
	int inited;

	MyFilterData() : noiseLevel(25), numPixels(4), pCleaner(NULL), markUnfiltered(false), inited(0) {}

	void init() {
		log("MyFilterData.init");
		noiseLevel = 25; numPixels = 4; markUnfiltered = false;
		inited = 1;
	}

	Cleaner *cln() {
		if (pCleaner==NULL) {
			pCleaner = new Cleaner(); 	
		}
		return pCleaner;
	}

	void deinit() {
		log("MyFilterData.deinit");
		if (pCleaner) {
			delete pCleaner;
			pCleaner = NULL;			
		}
	}

	void apply() {
		log("MyFilterData.apply");
		if (!pCleaner) {
			log("MyFilterData.apply: pCleaner is null! returning");
			return;
		}
		pCleaner->noiseAmplitude = noiseLevel;
		pCleaner->maxNoisyPixels = numPixels;
		pCleaner->markUnfiltered = markUnfiltered;
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

	return FILTERPARAM_SWAP_BUFFERS | FILTERPARAM_SUPPORTS_ALTFORMATS;// | FILTERPARAM_HAS_LAG(1);
}

sint64 prefetchProc(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, sint64 frame)
{
	return frame;
}

bool prefetchProc2(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, sint64 frame, IVDXVideoPrefetcher *prefetcher)
{
	int fr = frame;
	char str[256];
	sprintf(str, "prefetch fr=%d", fr);
	log(str);
	for(int f=fr-1; f<=fr+1; f++)
		prefetcher->PrefetchFrame(0, f, f);
	prefetcher->PrefetchFrameSymbolic(0, frame);
	return true;
}

int runProc(const VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("runProc");
	auto src = fa->src.mpPixmap;
	auto dst = fa->dst.mpPixmap;
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	Cleaner *c = pData->cln();
	if (c->inited==0) { 
		if (pData->inited==0)
			pData->init();
		pData->apply();
		c->init(src->w, src->h);
	}
	VDXFilterStateInfo *psi = fa->pfsi;
	int fn = 0;
	if (psi) {
		SHOW(psi->lCurrentFrame);
		fn = psi->lCurrentFrame;
	} else log("pfsi is null");
	SHOW(fa->mSourceFrameCount);
	assert(fa->mSourceFrameCount==3);
	for(int i=0;i<fa->mSourceFrameCount;i++)
		SHOW(fa->mpSourceFrames[i]->mFrameNumber);
	SHOW(fa->src.mFrameNumber);
	c->process(fa->mpSourceFrames, dst, fn);	
	return 0;
}

int startProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("startProc");
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	auto src = fa->src.mpPixmap;
	Cleaner *c = pData->cln();
	if (pData->inited==0)
		pData->init();
	pData->apply();
	c->init(src->w, src->h);
	return 0;
}

int endProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff) {
	log("endProc");
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	pData->deinit();
	return 0;
}

extern "C" __declspec(dllexport) void __cdecl benchmark(int w, int h, BYTE *frm0, BYTE *frm1, BYTE *frm2) {
	Cleaner c;
	c.init(w, h);
	VDXPixmap dst;
	dst.data = calloc(w*h,1);
	dst.data2 = calloc(w*h/4, 1);
	dst.data3 = calloc(w*h/4, 1);
	dst.pitch = w;
	dst.pitch2 = w/2;
	dst.pitch3 = w/3;
	dst.w = w;
	dst.h = h;

	BYTE *frms[3] = {frm0, frm1, frm2};
	VDXPixmap mps[3];
	VDXFBitmap bmps[3];
	for(int i=0;i<3;i++) {
		mps[i].pitch = w;
		mps[i].pitch2 = w/2;
		mps[i].pitch3 = w/2;
		mps[i].data = frms[i];
		mps[i].data2 = frms[i] + w*h;
		mps[i].data3 = frms[i] + w*h + w*h/4;
		mps[i].w = w;
		mps[i].h = h;
		bmps[i].mpPixmap = &mps[i];
	}

	VDXFBitmap *frames[3] = { &bmps[0], &bmps[1], &bmps[2] };
	for(int n=1;n<201;n++) {
		printf("\r%d..   ", n);
		int fr = n*5;
		bmps[0].mFrameNumber = fr - 1;
		bmps[1].mFrameNumber = fr;
		bmps[2].mFrameNumber = fr + 1;
		c.process(frames, &dst, fr);
	}

}

extern HINSTANCE g_hInst;

INT_PTR CALLBACK SettingsDlgProc(HWND hdlg, UINT msg, WPARAM wParam, LPARAM lParam) {
	MyFilterData* pData = (MyFilterData*)GetWindowLongPtr(hdlg, DWLP_USER);

    switch(msg) {
        case WM_INITDIALOG:
			log("init dialog");
            SetWindowLongPtr(hdlg, DWLP_USER, lParam);
			pData = (MyFilterData*)lParam;
			if (pData->inited==0) pData->init();

			SendMessage(GetDlgItem(hdlg, IDC_NOISELEVELSLIDER), TBM_SETRANGE, 0, MAKELONG(0, 100)); 
			SendMessage(GetDlgItem(hdlg, IDC_NOISELEVELSLIDER), TBM_SETPOS, 1, pData->noiseLevel); 
			SetDlgItemInt(hdlg, IDC_NOISELEVEL, pData->noiseLevel, FALSE);

			SendMessage(GetDlgItem(hdlg, IDC_COUNTSLIDER), TBM_SETRANGE, 0, MAKELONG(0, 60)); 
			SendMessage(GetDlgItem(hdlg, IDC_COUNTSLIDER), TBM_SETPOS, 1, pData->numPixels); 
			SetDlgItemInt(hdlg, IDC_COUNT, pData->numPixels, FALSE);

            CheckDlgButton(hdlg, IDC_MARKUNFILTERED, pData->markUnfiltered ? BST_CHECKED : BST_UNCHECKED);

			return TRUE;
		case WM_HSCROLL: {
			HWND h = (HWND)lParam;
			UINT id = GetWindowLong(h, GWL_ID);
			int x;
			switch(id) {
			case IDC_NOISELEVELSLIDER:
				x = SendMessage(h, TBM_GETPOS, 0, 0); 
				//pData->noiseLevel = x;
				SetDlgItemInt(hdlg, IDC_NOISELEVEL, x, FALSE);
				break;
			case IDC_COUNTSLIDER:
				x = SendMessage(h, TBM_GETPOS, 0, 0); 
				//pData->numPixels = x;
				SetDlgItemInt(hdlg, IDC_COUNT, x, FALSE);
				break;
			}
			SetWindowLong(hdlg, DWLP_MSGRESULT, 0);
			return TRUE;
		} break;
		case WM_COMMAND:
			switch(LOWORD(wParam)) {
				case IDOK:
					log("OK btn pressed");
					pData->noiseLevel = GetDlgItemInt(hdlg, IDC_NOISELEVEL, NULL, FALSE);
					pData->numPixels = GetDlgItemInt(hdlg, IDC_COUNT, NULL, FALSE);
					pData->markUnfiltered = IsDlgButtonChecked(hdlg, IDC_MARKUNFILTERED);
					pData->apply();
					EndDialog(hdlg, TRUE);
					return TRUE;
				case IDCANCEL:
					EndDialog(hdlg, FALSE);
					return TRUE;
				case IDC_HOMEPAGE:
					ShellExecute(NULL, NULL, L"http://www.infognition.com/", NULL, NULL, SW_SHOW);
					break;
			}
	}//switch msg
	return FALSE;
}

int configProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff, VDXHWND hwndParent) {
    MyFilterData* pData = (MyFilterData*)fa->filter_data;
	log("configProc");
	//pData->ifp = fa->ifp;
    auto res = !DialogBoxParam(g_hInst, MAKEINTRESOURCE(IDD_OPTIONS), (HWND)hwndParent, SettingsDlgProc, (LPARAM)pData);
    return res;
}

static void stringProc2(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, char *buf, int maxlen) {
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	_snprintf(buf, maxlen, " (%d, %d)", pData->noiseLevel, pData->numPixels );
}

void stringProc(const VDXFilterActivation *fa, const VDXFilterFunctions *ff, char *buf) {
	stringProc2(fa, ff, buf, 80);
}

void configScriptFunc(IVDXScriptInterpreter *isi, void *lpVoid, VDXScriptValue *argv, int argc) 
{
	VDXFilterActivation *fa = (VDXFilterActivation *)lpVoid;
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	if (pData->inited==0) pData->init();
	if (argc>=2) {
		pData->noiseLevel = argv[0].asInt(); 
		pData->numPixels = argv[1].asInt();
	}
}

VDXScriptFunctionDef script_functions[] = {
    { (VDXScriptFunctionPtr)configScriptFunc, "Config", "0ii" },
    { NULL, NULL, NULL },
};

VDXScriptObject script_obj = { NULL, script_functions, NULL };

bool fssProc(VDXFilterActivation *fa, const VDXFilterFunctions *ff, char *buf, int bufsize) 
{
	MyFilterData* pData = (MyFilterData*)fa->filter_data;
	_snprintf(buf, bufsize, "Config(%d, %d)", pData->noiseLevel, pData->numPixels);
	return true;
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
	&script_obj, fssProc,
	stringProc2,0,0,0,
	0, 0,
	prefetchProc2,0    
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
