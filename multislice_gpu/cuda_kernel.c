/*
    This version uses sequential addressing -- no divergence or bank conflicts.
*/
template <class T>
__global__ void reduce2(T *g_idata, T *g_odata, unsigned int n)
{
    T *sdata = SharedMemory<T>();

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

    sdata[tid] = (i < n) ? g_idata[i] : 0;

    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s=blockDim.x/2; s>0; s>>=1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }

        __syncthreads();
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

__global__ void module(float2 *g_idata, float *g_odata, uint32_t n) {

     unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
     if (i < n) {
        g_odata[i] = g_idata[i].x * g_idata[i].x + g_idata[i].y + g_idata[i].y;
     }

}
