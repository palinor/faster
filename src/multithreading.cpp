


#ifdef __APPLE__
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>

static pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
static const char *semaphore_name = "queue_semaphore";

typedef struct thread_container
{
	int JobId;
} thread_container;

typedef struct thread_queue
{
	thread_container *Containers;
	size_t MaxQueueSize;
	size_t volatile NextEntryToDo;
	size_t volatile EntryCompletionCount;
	size_t volatile EntryCount;
	sem_t *Semaphore;
} thread_queue;

void InitThreadQueue(thread_queue *Queue, size_t MaxQueueSize)
{
	Queue->Semaphore = sem_open(semaphore_name, O_CREAT);
	Queue->Containers = malloc(sizeof(thread_container) * MaxQueueSize);
	Queue->MaxQueueSize = MaxQueueSize;
	Queue->NextEntryToDo = 0;
	Queue->EntryCount = 0;
}

void PushContainerToQueue(thread_container *Container, thread_queue *Queue, sem_t *Semaphore)
{
	memcpy(Queue->Containers + Queue->EntryCount++, Container, sizeof(thread_container));
	sem_post(Semaphore);
}

typedef struct thread_input
{
	thread_queue *Queue;
	int ThreadIdx;
} thread_input;

bool IsWorkRemaining(thread_queue *Queue)
{
	return (Queue->EntryCompletionCount != Queue->EntryCount);
}

typedef struct work_queue_entry
{
	size_t EntryIdx;
	bool IsValid;
} work_queue_entry;

work_queue_entry GetNextElement(thread_queue *Queue)
{
	work_queue_entry Result;
	Result.IsValid = 0;
	if (Queue->NextEntryToDo < Queue->EntryCount)
	{
		pthread_mutex_lock(&queue_mutex);
		Result.EntryIdx = Queue->NextEntryToDo++;
		pthread_mutex_unlock(&queue_mutex);
		Result.IsValid = 1;
	}
	return Result;
}

void MarkEntryComplete(thread_queue *Queue, work_queue_entry Entry)
{
	pthread_mutex_lock(&queue_mutex);
	Queue->EntryCompletionCount++;
	pthread_mutex_unlock(&queue_mutex);
}

bool DoThreadWork(thread_queue *Queue, size_t thread_idx)
{
	work_queue_entry Entry = GetNextElement(Queue);
	if (Entry.IsValid)
	{
		thread_container *Container = Queue->Containers + Entry.EntryIdx;
		printf("Wrote job %d from thread %d\n", (int)Container->JobId, (int)thread_idx);
		MarkEntryComplete(Queue, Entry);
	}
	return Entry.IsValid;
}

void HandleThread(void *Input)
{
	thread_input *InputArgs = (thread_input *)Input;
	for (;;)
	{
		if (!DoThreadWork(InputArgs->Queue, InputArgs->ThreadIdx))
		{
			sem_wait(InputArgs->Queue->Semaphore);
		}
	}
}

int main()
{
	size_t NThreads = 8;
	size_t NJobs = 32;
	thread_queue Queue;
	InitThreadQueue(&Queue, 256);
	thread_input ThreadInput[NThreads];
	pthread_t MyThreads[NThreads];
	sem_t *Semaphore = sem_open(semaphore_name, O_CREAT);
	for (size_t ThreadIdx = 0; ThreadIdx < NThreads; ThreadIdx++)
	{
		ThreadInput[ThreadIdx].Queue = &Queue;
		ThreadInput[ThreadIdx].ThreadIdx = ThreadIdx;
		pthread_create(MyThreads + ThreadIdx, NULL, (void *)&HandleThread, (void *)(ThreadInput + ThreadIdx));
	}
	thread_container Container;
	for (size_t JobId = 0; JobId < NJobs; JobId++)
	{
		Container.JobId = JobId;
		PushContainerToQueue(&Container, &Queue, Semaphore);
	}
	while (IsWorkRemaining(&Queue))
	{
		DoThreadWork(&Queue, NThreads);
	}

	return 0;
}
#endif
