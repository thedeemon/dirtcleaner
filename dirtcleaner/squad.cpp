#include "stdafx.h"
#include "squad.h"
#include <malloc.h>

DWORD WINAPI SquadWorkerThreadProc(LPVOID lpParameter)
{
	((CSquadWorker*)lpParameter)->ThreadProc();
	return 0;
}

void CSquadWorker::GetSegment(int totalsize, int &segstart, int &segsize)
{
	if (totalsize >= pSquad->nw) {
		segstart = totalsize * myNum / pSquad->nw;
		int segend = totalsize * (myNum+1) / pSquad->nw;
		if (segend > totalsize)
			segend = totalsize;
		segsize = segend - segstart;
	} else { //totalsize < nw
		if (myNum < totalsize) {
			segstart = myNum; segsize = 1;
		} else {
			segstart = segsize = 0;
		}
	}
}

/////////////////////////////////////////////////////////////////////

CSquad::CSquad(int nThreads)
{
	if (nThreads<1) 
		nThreads = 1;
	nw = nThreads;
	workers = (CSquadWorker**)calloc(nThreads, sizeof(CSquadWorker*));
	for(int i=0;i<nThreads; i++)
		workers[i] = new CSquadWorker(this, i);
	if (nThreads>1) {
		ev_free = (HANDLE*)calloc(nThreads, sizeof(HANDLE));
		ev_havejob = (HANDLE*)calloc(nThreads, sizeof(HANDLE));
		ev_sync = (HANDLE*)calloc(nThreads, sizeof(HANDLE));
		DWORD tid = 0;
		for(int i=0;i<nThreads; i++) {
			ev_free[i] =  CreateEvent(NULL, TRUE/*manual*/, FALSE/*initial*/, NULL);
			ev_havejob[i] = CreateEvent(NULL, FALSE/*auto*/, FALSE, NULL);
			ev_sync[i] =    CreateEvent(NULL, FALSE/*auto*/, FALSE, NULL);
			workers[i]->thread_handle = CreateThread(NULL,256*1024, SquadWorkerThreadProc, workers[i], 0, &tid);
		}
	}
}

CSquad::~CSquad()
{
	if (nw>1) { //stop threads
		WaitTillAllFree();
		cur_command = -1;		
		SignalJob();
		HANDLE* ths = (HANDLE*)calloc(nw, sizeof(HANDLE));
		for(int i=0;i<nw;i++)
			ths[i] = workers[i]->thread_handle;
		WaitForMultipleObjects(nw, ths, TRUE, INFINITE);
		free(ths);

		for(int i=0;i<nw;i++) {
			CloseHandle(ev_free[i]);
			CloseHandle(ev_havejob[i]);
			CloseHandle(ev_sync[i]);
		}
		free(ev_free); free(ev_havejob); free(ev_sync);
	}
	for(int i=0;i<nw; i++) {
		if (workers[i]->thread_handle)
			CloseHandle(workers[i]->thread_handle);
		delete workers[i];
	}
	free(workers);
}

void CSquad::WaitTillAllFree()
{
	if (nw<2) return;
	WaitForMultipleObjects(nw, ev_free, TRUE/*all*/, INFINITE);
}

void CSquad::SignalJob() //signal to workers they have a job
{
	if (nw<2) return;
	for(int i=0;i<nw;i++)
		SetEvent(ev_havejob[i]);
}

void CSquad::Sync(int mynum)
{
	if (nw<2) return;	
	if (mynum > 0) {
		SignalObjectAndWait(ev_sync[mynum], ev_havejob[mynum], INFINITE, FALSE);
	} else {
		WaitForMultipleObjects(nw-1, &ev_sync[1], TRUE/*all*/, INFINITE);
		for(int i=1;i<nw;i++)
			SetEvent(ev_havejob[i]);
	}
}

void CSquad::RunParallel(int command, void *params, ISquadJob *job)
{
	cur_command = command;
	cur_params = params;
	cur_job = job;

	if (nw>1) { //run in worker threads
		WaitTillAllFree();
		for(int i=0;i<nw;i++)
			ResetEvent(ev_free[i]);
		SignalJob();
		WaitTillAllFree();
	} else  //run in this thread
		job->RunCommand(command, params, workers[0]);
}

void CSquad::ThreadProc(CSquadWorker *sqworker)
{	
	const int mynum = sqworker->MyNum();
	while(1) {
		DWORD waitres = SignalObjectAndWait(ev_free[mynum], ev_havejob[mynum], INFINITE, FALSE);
		if (waitres==WAIT_OBJECT_0) {//signaled
			if (cur_command < 0) {//stop
				break; 
			}
			cur_job->RunCommand(cur_command, cur_params, sqworker);
		} 
	}
}
