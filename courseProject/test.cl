kernel void MultiplicationMatrix(global float* mA,
								 global float* mB,
								 global float* mC, 
								 private unsigned int sz)
{
	int row = get_global_id(0);
	int col = get_global_id(1);
	float temp = 0;
	for (int i = 0; i < sz; ++i) {
		temp += mA[row * sz + i] * mB[i * sz + col];
	}
	mC[row * sz + col] = temp;
}

