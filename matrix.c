#include <assert.h>
#include <pthread.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

typedef struct Matrixf64
{
	size_t NRows;
	size_t NCols;
	double *Data;
} Matrixf64;

/*
Originally templated off of Savine's intro to Modern Computational Finance.
*/
Matrixf64 Matrixf64MultiplyFast(const Matrixf64 *A, const Matrixf64 *B)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	for (size_t i = 0; i < A->NRows; i++)
	{
		const double *ARow = A->Data + i * A->NCols;
		double *ResultRow = Result.Data + i * Result.NCols;
		for (size_t k = 0; k < B->NRows; k++)
		{
			const double *BRow = B->Data + k * B->NCols;
			const double Aik = ARow[k];
			// #pragma clang loop vectorize_width(2) interleave_count(8)
			for (size_t j = 0; j < B->NCols; j++)
			{
				ResultRow[j] += Aik * BRow[j];
			}
		}
	}
	return Result;
}

/*
Naive implementation of matrix multiplication. The thing is, there does not seem to be a noticeable difference
between this and the fast one above when turning on the -Ofast flag in clang.

The difference is noticeable at lower levels of optimization though. Some quick benchmarks for multiplying 2 random 1000 x 1000
matrices:

-O0:
Slow done in 3.454435s
Fast done in 1.221509s

-O1:
Slow done in 3.452416s
Fast done in 0.326448s

-O2:
Slow done in 0.841188s
Fast done in 0.158153s

-O3:
Slow done in 0.843141s
Fast done in 0.159574s

-Ofast:
Slow done in 0.163217s
Fast done in 0.156944s

Moral of the story, compilers are good at optimizing shit.

Numpy benchmark: running python 3.10, numpy 1.22.3, multiplying 2 1000x1000 random matrices is done in 0.0199429988861084s....

So python is still better. We need multithreading! With compiler optimizations on, we are already better than the numbers Savine reports
in his book as benchmarks on his Intel (on a MacBook Pro)

Pytorch 1.11.0:
Done in 0.0027608871459960938s

An order of magnitude faster than numpy?

*/

Matrixf64 Matrixf64MultiplyNaive(const Matrixf64 *A, const Matrixf64 *B)
{
	assert(A->NCols == B->NRows);
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	for (size_t i = 0; i < A->NRows; i++)
	{
		for (size_t j = 0; j < B->NCols; j++)
		{
			for (size_t k = 0; k < B->NRows; k++)
			{
				Result.Data[i * A->NRows + j] += A->Data[i * A->NRows + k] * B->Data[j * B->NRows + k];
			}
		}
	}
	return Result;
}

/*
Let's take the fast version and make it multithreaded

*/

typedef struct matrix_multiply_job
{
	double *ARow;
	double *BRow;
	size_t ARowIdx;
	size_t BRowIdx;
	size_t ARowLen;
	size_t BRowLen;
	volatile double *ResultRow;
	pthread_mutex_t *ResultRowLock;
} matrix_multiply_job;

typedef struct thread_queue
{
	size_t MaxQueueSize;
	size_t volatile NextEntryToDo;
	size_t volatile EntryCompletionCount;
	size_t volatile EntryCount;
	sem_t *Semaphore;
	pthread_mutex_t Lock;
} thread_queue;

void Matrixf64Multiply2Rows(matrix_multiply_job *Entry)
{
	pthread_mutex_lock(Entry->ResultRowLock);
	#pragma clang loop vectorize(enable)
	for (size_t k = 0; k < Entry->BRowLen; ++k)
	{
		// we just need to lock this so only one thread can write to it at a time
		Entry->ResultRow[k] += Entry->ARow[k] * Entry->BRow[k];
	}
	pthread_mutex_unlock(Entry->ResultRowLock);
}

void InitThreadQueue(thread_queue *Queue, size_t MaxQueueSize)
{
	pthread_mutex_init(&(Queue->Lock), NULL);
	Queue->Semaphore = sem_open("queue_semaphore", O_CREAT);
	Queue->MaxQueueSize = MaxQueueSize;
	Queue->NextEntryToDo = 0;
	Queue->EntryCount = 0;
	Queue->EntryCompletionCount = 0;
}

void PushContainerToQueue(thread_queue *Queue)
{
	Queue->EntryCount++;
	sem_post(Queue->Semaphore);
}

typedef struct thread_input
{
	thread_queue *Queue;
	matrix_multiply_job *Jobs;
	int ThreadIdx;
} thread_input;

bool IsWorkRemaining(thread_queue *Queue)
{
	pthread_mutex_lock(&(Queue->Lock));
	bool WorkLeft = Queue->EntryCompletionCount < Queue->EntryCount;
	pthread_mutex_unlock(&(Queue->Lock));
	return WorkLeft;
}

typedef struct work_queue_entry
{
	size_t EntryIdx;
	bool IsValid;
} work_queue_entry;

work_queue_entry GetNextElement(thread_queue *Queue)
{
	work_queue_entry Result;
	Result.IsValid = false;
	pthread_mutex_lock(&(Queue->Lock));
	if (Queue->NextEntryToDo < Queue->EntryCount)
	{
		Result.EntryIdx = Queue->NextEntryToDo++;
		Result.IsValid = true;
	}
	pthread_mutex_unlock(&(Queue->Lock));
	return Result;
}

void MarkEntryComplete(thread_queue *Queue)
{
	pthread_mutex_lock(&(Queue->Lock));
	Queue->EntryCompletionCount++;
	pthread_mutex_unlock(&(Queue->Lock));
}

bool DoMatrixMultiplyf64ThreadWork(matrix_multiply_job *Jobs, thread_queue *Queue)
{
	work_queue_entry Entry = GetNextElement(Queue);
	if (Entry.IsValid)
	{
		matrix_multiply_job *Job = Jobs + Entry.EntryIdx;
		Matrixf64Multiply2Rows(Job);
		MarkEntryComplete(Queue);
	}
	return Entry.IsValid;
}

void HandleMatrixMultiplyThread(void *Input)
{
	thread_input *InputArgs = (thread_input *)Input;
	for (;;)
	{
		if (IsWorkRemaining(InputArgs->Queue)) {
			DoMatrixMultiplyf64ThreadWork(InputArgs->Jobs, InputArgs->Queue);
		} else{
			sem_wait(InputArgs->Queue->Semaphore);
		}
	}
}


Matrixf64
Matrixf64MultiplyMultithreaded(const Matrixf64 *A, const Matrixf64 *B, size_t NThreads)
{
	assert(A->NCols == B->NRows);
	thread_input ThreadInfo[NThreads];
	matrix_multiply_job *Jobs = malloc(sizeof(matrix_multiply_job) * A->NRows * B->NRows);
	pthread_t Threads[NThreads];
	thread_queue *Queue = malloc(sizeof(thread_queue));
	InitThreadQueue(Queue, 1000 * 1000 + 1);
	for (size_t i = 0; i < NThreads; i++)
	{
		ThreadInfo[i].ThreadIdx = i;
		ThreadInfo[i].Jobs = Jobs;
		ThreadInfo[i].Queue = Queue;
		pthread_create(Threads + i, NULL, (void *)&HandleMatrixMultiplyThread, (void *)(ThreadInfo + i));
	}
	Matrixf64 Result;
	Result.NRows = A->NRows;
	Result.NCols = B->NCols;
	Result.Data = calloc(Result.NRows * Result.NCols, sizeof(double));
	assert(A->NRows * B->NRows < Queue->MaxQueueSize);
	pthread_mutex_t RowLocks[A->NRows];
	for (size_t i = 0; i < A->NRows; i++){
		pthread_mutex_init(RowLocks + i, NULL);
	}
	for (size_t j = 0; j < B->NRows; j++)
	{
		for (size_t i = 0; i < A->NRows; i++)
		{	
			// go in this order so that we stack the jobs in column major of the result matrix -> maybe better order for the multithreading queue
			matrix_multiply_job *ThisJob = Jobs + j * B->NRows + i;
			ThisJob->ARowIdx = i;
			ThisJob->BRowIdx = j;
			ThisJob->ARow = A->Data + i * A->NCols;
			ThisJob->BRow = B->Data + j * B->NCols;
			ThisJob->ARowLen = A->NCols;
			ThisJob->BRowLen = B->NCols;
			ThisJob->ResultRow = Result.Data + i * A->NCols;
			ThisJob->ResultRowLock = RowLocks + i;
			PushContainerToQueue(Queue);
		}
	}
	while (IsWorkRemaining(Queue))
	{
		DoMatrixMultiplyf64ThreadWork(Jobs, Queue);
	}
	InitThreadQueue(Queue, 0); // set queue size back to 0 so all threads go back to sleep when funciton returns
	free(Jobs);
	return Result;
}

Matrixf64 RandomMatrixf64(size_t NRows, size_t NCols)
{
	srand(0);
	Matrixf64 Result;
	Result.NRows = NRows;
	Result.NCols = NCols;
	Result.Data = malloc(NRows * NCols * sizeof(double));
	for (size_t i = 0; i < NRows; i++)
	{
		for (size_t j = 0; j < NCols; j++)
		{
			Result.Data[i * NRows + j] = ((double)rand()) / RAND_MAX;
		}
	}
	return Result;
}

Matrixf64 Identity(size_t MatrixSize)
{
	Matrixf64 Result;
	Result.NRows = MatrixSize;
	Result.NCols = MatrixSize;
	Result.Data = calloc(MatrixSize * MatrixSize, sizeof(double));
	for (size_t i = 0; i < MatrixSize; i++)
	{
		Result.Data[i * MatrixSize + i] = 1;
	}
	return Result;
}

Matrixf64 Ones(size_t MatrixRows, size_t MatrixCols)
{
	Matrixf64 Result;
	Result.NRows = MatrixRows;
	Result.NCols = MatrixCols;
	Result.Data = malloc(MatrixRows * MatrixCols * sizeof(double));
	for (size_t i = 0; i < MatrixRows; i++)
	{
		for (size_t j = 0; j < MatrixCols; j++)
		{
			Result.Data[i * MatrixCols + j] = 1;
		}
	}
	return Result;
}

/* This is technically the slow version, but we get a x20 speedup with clang on -Ofast */
void test_matrix_multiply(size_t MatrixSize, int PrintResult, int NThreads, bool UseRandom)
{
	Matrixf64 A, B;
	if (UseRandom)
	{
		A = RandomMatrixf64(MatrixSize, MatrixSize);
		B = RandomMatrixf64(MatrixSize, MatrixSize);
	}
	else
	{
		A = Ones(MatrixSize, MatrixSize), B = Identity(MatrixSize);
	}
	clock_t start = clock();
	Matrixf64 Result = Matrixf64MultiplyNaive(&A, &B);
	clock_t end = clock();
	free(Result.Data);
	printf("Slow done in %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);
	start = clock();
	Result = Matrixf64MultiplyFast(&A, &B);
	end = clock();
	printf("Fast done in %lfs\n", (double)(end - start) / CLOCKS_PER_SEC);
	free(Result.Data);
	start = clock();
	Result = Matrixf64MultiplyMultithreaded(&A, &B, NThreads);
	end = clock();
	printf("Multithreaded done in %lfs with %d threads\n", (double)(end - start) / CLOCKS_PER_SEC, NThreads);
	if (PrintResult)
	{
		printf("Result\n");
		for (size_t i = 0; i < Result.NRows; i++)
		{
			for (size_t j = 0; j < Result.NCols; j++)
			{
				printf("%f ", Result.Data[i * Result.NRows + j]);
			}
			printf("\n");
		}
	}
}

int main()
{
	test_matrix_multiply(1000, 0, 8, 0);
	return 0;
}
