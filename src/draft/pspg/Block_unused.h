

   static __global__ 
   void pointwiseMul(const cudaReal* a, 
                     const cudaReal* b, 
                     cudaReal* result, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < size; i += nThreads) {
         result[i] = a[i] * b[i];
      }
   }

   static __global__ 
   void pointwiseFloatMul(const cudaReal* a, 
                          const float* b, 
                          cudaReal* result, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < size; i += nThreads) {
         result[i] = a[i] * b[i];
         // printf("result[%d], =  %d\n", i , result[i]);
      }
   }

   static __global__ 
   void equalize (const cudaReal* a, double* result, int size)
   {  
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < size; i += nThreads) {
         result [i] = a [i];
      }
   }

   static __global__ 
   void pointwiseMulUnroll2(const cudaReal* a, const cudaReal* b, 
                            cudaReal* result, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x * 2 + threadIdx.x * 2;
      cudaReal localResult[2];
      for (int i = startID; i < size; i += nThreads * 2) {
         localResult[0] = a[i] * b[i];
         localResult[1] = a[i + 1] * b[i + 1];
         result[i] = localResult[0];
         result[i + 1] = localResult[1];
         //result[i] = a[i] * b[i];
         //result[i + 1] = a[i + 1] * b[i + 1];

      }
   }

   static __global__ 
   void pointwiseMulCombi(cudaReal* a,const cudaReal* b, 
                          cudaReal* c,const cudaReal* d,
                          const cudaReal* e, int size) 
   {
      //c = a * b
      //a = d * e
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      cudaReal tempA;
      for (int i = startID; i < size; i += nThreads) {
         tempA = a[i];
         c[i] = tempA * b[i];
         a[i] = d[i] * e[i];

      }
   }

   static __global__ 
   void richardsonExp(cudaReal* qNew, const cudaReal* q1, 
                      const cudaReal* q2, int size) 
   {
      int nThreads = blockDim.x * gridDim.x;
      int startID = blockIdx.x * blockDim.x + threadIdx.x;
      for (int i = startID; i < size; i += nThreads) {
         qNew[i] = (4.0 * q2[i] - q1[i]) / 3.0;
      }
   }

