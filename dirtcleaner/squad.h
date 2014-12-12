#ifndef _SQUAD_H_
#define _SQUAD_H_

//#if WINVER < 0x0400
#define WINVER 0x0400
#define _WIN32_WINNT 0x0400
//#endif

#include <windows.h>

class CSquadWorker;


class ISquadJob {
public:
	virtual void RunCommand(int command, void *params, CSquadWorker *sqworker)=0 ;
};

class CSquad {
	friend class CSquadWorker;

	CSquadWorker **workers;
	int nw; //number of workers
	//long active_workers, unsync_workers;
	//HANDLE job_event, jobdone_event, sync_event, sync_mutex, free_sem;
	HANDLE* ev_free;
	HANDLE* ev_havejob;
	HANDLE* ev_sync;

	int cur_command;
	void *cur_params;
	ISquadJob *cur_job;
	

	void Sync(int mynum);
	void ThreadProc(CSquadWorker *sqworker);
	void WaitTillAllFree();
	void SignalJob(); //signal to workers they have a job

public:
	CSquad(int nThreads);
	~CSquad();

	int NumThreads() { return nw; }
	void RunParallel(int command, void *params, ISquadJob *job);
};

class CSquadWorker {
	CSquad *pSquad;
	int myNum;

public:
	HANDLE thread_handle;

	//called from Squad
	CSquadWorker(CSquad *squad, int mynum) : pSquad(squad), myNum(mynum), thread_handle(NULL) {}
	void ThreadProc() {	pSquad->ThreadProc(this); };

	//called from Job
	void GetSegment(int totalsize, int &segstart, int &segsize);
	void Sync()  { pSquad->Sync(myNum); }
	int NumThreads() { return pSquad->NumThreads(); }
	int MyNum() { return myNum; }

};


#endif