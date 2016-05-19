#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

kernel void MultMatrix(global float* mA,
                       global float* mB,
                       global float* mC,
                       private unsigned int sz)
{
    int row = get_global_id(0);
    int col = row - row / sz * sz;
    row /= sz;
    float temp = 0;
    for (int k = 0; k < sz; ++k)
    {
        temp += mA[row * sz + k] * mB[k * sz + col];
    }
    mC[row * sz + col] = temp;
}

kernel void FindMaxABS(global float* a, 
					   global float* maxElem, 
					   global int* posMax,
					   private unsigned int sz)
{
    int pos = get_global_id(0);
    int col = 0;
    int row = pos * sz;
    float _max = -1;
	int posI, posJ;
    for (int i = 0; i < sz; ++i) {
        if (i == pos) continue;
        if (_max < fabs(a[row + i])) {
            _max = fabs(a[row + i]);
			posJ = i;
        }
    }
    maxElem[pos] = _max;
    posMax[pos] = posJ;
    barrier(CLK_GLOBAL_MEM_FENCE);
    if (!pos) {
    	_max = -1;
    	posJ = -1;
    	for (int i = 0; i < sz; ++i) {
    		if (_max < maxElem[i]) {
    			_max = maxElem[i];
    			posI = i;
    			posJ = posMax[i];
    		}
    	}
    	posMax[0] = posI;
    	posMax[1] = posJ;
    }
}

kernel void ComputeNorm(global float* a,
                        global float* res,
			      		private unsigned int sz)
{
    int pos = get_global_id(0);
    int row = pos * sz;
    float sum = 0;
    for (int i = 0; i < sz; ++i) {
        if (pos == i) continue;
        sum += a[row + i] * a[row + i];
    }
    res[pos] = sum;
    barrier(CLK_GLOBAL_MEM_FENCE);
    if (!pos) {
    	sum = 0;
    	for (int i = 0; i < sz; ++i) {
    		sum += res[i];
    	}
    	res[0] = sum;
    }
}

kernel void TransposeMatrix(global float* a,
                            global float* res,
                            private unsigned int sz)
{
    int pos = get_global_id(0);
    int row = pos * sz;
    for (int i = 0; i < sz; ++i) {
        res[pos + i * sz] = a[row + i]; 
    }
}
