#include <BinSearch.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif
#include <common.h>

using namespace BinSearch;

// Cross-platform thread handling

#ifdef _WIN32
typedef HANDLE thread_t;

typedef struct ThreadData {
    void* (*start_routine)(void*);
    void* arg;
} ThreadData;

DWORD WINAPI ThreadProc(LPVOID lpParam) {
    ThreadData* data = (ThreadData*)lpParam;
    data->start_routine(data->arg);
    free(data); // Clean up after using
    return 0;
}

void create_thread(thread_t* thread, void* (*start_routine) (void*), void* arg) {
    ThreadData* data = (ThreadData*)malloc(sizeof(ThreadData));
    data->start_routine = start_routine;
    data->arg = arg;
    *thread = CreateThread(NULL, 0, ThreadProc, data, 0, NULL);
}

void join_thread(thread_t thread) {
    WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
}

#else
typedef pthread_t thread_t;

void create_thread(thread_t* thread, void* (*start_routine) (void*), void* arg) {
    pthread_create(thread, NULL, start_routine, arg);
}

void join_thread(thread_t thread) {
    pthread_join(thread, NULL);
}
#endif

void dequantize_cpu(float *code, unsigned char *A, float *absmax, float *out, long long blocksize, long long n) {
    for (long long block_idx = 0; block_idx < n; block_idx += blocksize) {
        long long valid_items = n - block_idx >= blocksize ? blocksize : n - block_idx;
        long long block_end = block_idx + valid_items;
        for (long long i = block_idx; i < block_end; i++)
            out[i] = code[A[i]] * absmax[block_idx / blocksize];
    }
}

void quantize_cpu(float *code, float *A, float *absmax, unsigned char *out, long long blocksize, long long n)
{

    // the default code is has range [-0.993, 1.0] which can cause an error in the binary search algorithm used below
    code[0] = -1.0f;

    long long num_blocks = n / blocksize;
    num_blocks += n % blocksize == 0 ? 0 : 1;

    const uint32 elements_code = 256;
    BinAlgo<Scalar, float, Direct2> bin_searcher(code, elements_code);

    int thread_wave_size = 256;
    // we chunk the thresds into waves of 256 since the max limit is
    // between 16k and 64k on Linux (we reach this when running BLOOM-176B with a large batch size)
    for(long long offset = 0; offset < num_blocks; offset+=thread_wave_size)
    {
      long long valid_chunks = num_blocks - offset >= thread_wave_size ? thread_wave_size : num_blocks - offset;
      thread_t* threads = (thread_t*)malloc(sizeof(thread_t) * valid_chunks);

      struct quantize_block_args **args = (quantize_block_args **) malloc(valid_chunks * sizeof(quantize_block_args *));

      for(long long i = 0; i < valid_chunks; i++)
          args[i] = (quantize_block_args *) malloc(sizeof(quantize_block_args));

      int chunks_processed = 0;
      for(long long block_idx = offset*blocksize; block_idx < n; block_idx += blocksize)
      {
          long long valid_items = n - block_idx >= blocksize ? blocksize : n - block_idx;
          long long block_end = block_idx + valid_items;

          struct quantize_block_args *arg = args[chunks_processed];
          arg->bin_searcher = &bin_searcher;
          arg->code = code;
          arg->A = A;
          arg->absmax = absmax;
          arg->out = out;
          arg->block_end = block_end;
          arg->block_idx = block_idx;
          arg->threadidx = block_idx / blocksize;
          arg->blocksize = blocksize;

          create_thread(&threads[chunks_processed], &quantize_block, (void*)arg);
          chunks_processed += 1;
          if(chunks_processed == valid_chunks){ break; }
      }

      for (int i = 0; i < valid_chunks; i++)
          join_thread(threads[i]);

      free(threads);
      for (int i = 0; i < valid_chunks; i++)
          free(args[i]);
      free(args);

    }

}
