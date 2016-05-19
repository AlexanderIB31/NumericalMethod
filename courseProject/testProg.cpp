#include <CL/cl.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <chrono>
#include <cstring>

#define MAX_SOURCE_SIZE (0x100000)

using namespace std;

void checkRet(const cl_int err, const int raw) {
    if (err != CL_SUCCESS) {
        std::cerr << "Trouble: " << err << " | " << raw << std::endl;
        exit(1);
    }
}

cl_float Off(cl_float* matrix, cl_uint szMatrix, cl_kernel kernel,
             cl_command_queue command_queue, cl_context context, size_t* global_work_size) {
    cl_int err;
    szMatrix = max(1u, szMatrix);
    cl_float* norms     = (cl_float*)malloc(sizeof(cl_float) * szMatrix);
    cl_mem buffMatrix = clCreateBuffer(context, CL_MEM_READ_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
    checkRet(err, __LINE__);
    cl_mem buffNorm = clCreateBuffer(context, CL_MEM_READ_WRITE, szMatrix * sizeof(cl_float), NULL, &err);
    checkRet(err, __LINE__);    
    err = clEnqueueWriteBuffer(command_queue, buffMatrix, CL_FALSE, 0, szMatrix * szMatrix * sizeof(cl_float), matrix, 0, NULL, NULL);        
    checkRet(err, __LINE__);
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buffMatrix);
    checkRet(err, __LINE__);
    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&buffNorm);
    checkRet(err, __LINE__);
    err = clSetKernelArg(kernel, 2, sizeof(cl_uint), (void *)&szMatrix);
    checkRet(err, __LINE__);    
    err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
    checkRet(err, __LINE__);
    err = clEnqueueReadBuffer(command_queue, buffNorm, CL_TRUE, 0, szMatrix * sizeof(cl_float), norms, 0, NULL, NULL);
    checkRet(err, __LINE__);
    clReleaseMemObject(buffMatrix);    
    clReleaseMemObject(buffNorm);
    cl_float norm = norms[0];
    free(norms); 
    return sqrt(norm);   
}

void RotateMethodGPU() {
    string plotFileGPU = "plotDataGPU";
    // инициализация устройств
    // ************************************************************************************
    cl_int              err;
    cl_platform_id*     platforms;
    cl_uint             num_platforms;
    cl_device_id*       devices;
    cl_uint             num_devices;
    cl_context          context;
    cl_command_queue    command_queue;
    cl_program          program                     = NULL;
    cl_kernel           kernelFindMaxAbs            = NULL;
    cl_kernel           kernelMultMatrix            = NULL;
    cl_kernel           kernelComputeNorm           = NULL;
    cl_kernel           kernelTranspose             = NULL;
    cl_context_properties platforms_properties[3];

    /* получить доступные платформы */
    err = clGetPlatformIDs(0, NULL, &num_platforms);
    checkRet(err, __LINE__);
    platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id) * num_platforms);
    err = clGetPlatformIDs(num_platforms, platforms, NULL);
    checkRet(err, __LINE__);
    /* получить доступные устройства */
    // цикл нужен если платформ с GPU > 1
    //for (size_t i = 0; i < num_platforms; ++i) {
        platforms_properties[0] = (cl_context_properties)CL_CONTEXT_PLATFORM;  // indicates that next element is platform
        platforms_properties[1] = (cl_context_properties)platforms[0];  // platform is of type cl_platform_id
        platforms_properties[2] = (cl_context_properties)0;   // last element must b
        err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
        checkRet(err, __LINE__);
        devices = (cl_device_id*)malloc(sizeof(cl_device_id) * num_devices);
        err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
        checkRet(err, __LINE__);
    //}
    size_t param_sz;
    cl_uint param_int, valFPS = 0;
    char* param_string;
    for (size_t i = 0; i < num_devices; ++i) {        
        err = clGetDeviceInfo(devices[i], CL_DEVICE_ADDRESS_BITS, sizeof(cl_uint), (void*)&param_int, NULL);
        checkRet(err, __LINE__);        
        cout << "ADDRESS BITS = " << param_int << endl;
        err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), (void*)&param_int, NULL);
        checkRet(err, __LINE__);        
        valFPS = param_int;
        err = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint), (void*)&param_int, NULL);
        checkRet(err, __LINE__);        
        valFPS *= param_int;
        err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 0, NULL, &param_sz);
        checkRet(err, __LINE__);        
        param_string = (char*)malloc(sizeof(char) * param_sz);
        err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, param_sz, (void*)param_string, NULL);
        checkRet(err, __LINE__);
        cout << "FPS = " << valFPS << " | " << "device name: " << param_string << endl;
        free(param_string);
        err = clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, 0, NULL, &param_sz);
        checkRet(err, __LINE__);        
        param_string = (char*)malloc(sizeof(char) * param_sz);
        err = clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, param_sz, (void*)param_string, NULL);
        checkRet(err, __LINE__);
        cout << "vendor: " << param_string << endl;
        cout << "**********************************" << endl;
        free(param_string);
    }
    // 3-ий аргумен = 0, т.к. всего одно устройство
    context = clCreateContextFromType(platforms_properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &err);    
    checkRet(err, __LINE__);
    command_queue = clCreateCommandQueue(context, devices[0], 0, &err);
    checkRet(err, __LINE__);

    FILE* fp;
    const char* fileName = "./test.cl";
    char *source_str;
    size_t source_size;

    try {
        fp = fopen(fileName, "r");
        if (!fp) {
            std::cerr << "Failed to load kernelMultMatrix." << std::endl;
            exit(1);
        }

        source_str = (char *)malloc(MAX_SOURCE_SIZE);
        source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
        fclose(fp);
    }
    catch (int a) {
        printf ("%d", a);
    }
    /* создать бинарник из кода программы */
    program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &err);
    checkRet(err, __LINE__);    
    /* скомпилировать программу */    
    err = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
    // build failed
    if (err != CL_SUCCESS) {
        cl_build_status status;
        size_t logSize;
        // check build error and build status first
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_STATUS, 
                sizeof(cl_build_status), &status, NULL);
 
        // check build log
        clGetProgramBuildInfo(program, devices[0], 
                CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        char* programLog = (char*)malloc((logSize + 1) * sizeof(char));
        clGetProgramBuildInfo(program, devices[0], 
                CL_PROGRAM_BUILD_LOG, logSize + 1, programLog, NULL);
        printf("Build failed; error=%d, status=%d, programLog:nn%s", 
                err, status, programLog);
        free(programLog);
    }
    // ************************************************************************************
    cl_float EPS = 0.01;
    cl_float EPSk;
    
    /*  создать кернел */
    kernelMultMatrix = clCreateKernel(program, "MultMatrix", &err);
    checkRet(err, __LINE__);
    kernelFindMaxAbs = clCreateKernel(program, "FindMaxABS", &err);
    checkRet(err, __LINE__);
    kernelMultMatrix = clCreateKernel(program, "MultMatrix", &err);
    checkRet(err, __LINE__);
    kernelComputeNorm = clCreateKernel(program, "ComputeNorm", &err);
    checkRet(err, __LINE__);
    kernelTranspose = clCreateKernel(program, "TransposeMatrix", &err);
    checkRet(err, __LINE__);

    cl_float*   matrixA         = NULL;
    cl_float*   matrixB         = NULL;
    cl_float*   matrixC         = NULL;
    cl_float*   matrixRotate    = NULL;
    cl_float*   matrixRotateT   = NULL;
    
    ofstream out(plotFileGPU, ios::out);
    for (cl_uint cnt = 3; cnt <= 3; cnt += 50) {
        cl_uint szMatrix = cnt;
        size_t global_work_size[2] = { szMatrix, szMatrix };

        // квадратная матрица, симметричная
        matrixA         = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        matrixB         = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        matrixC         = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        matrixRotate    = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        matrixRotateT   = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);        

        cl_mem      buffMatrixA      = NULL;
        cl_mem      buffMatrixB      = NULL;
        cl_mem      buffMatrixC      = NULL;

        for (cl_uint i = 0; i < szMatrix; ++i) {
            for (cl_uint j = 0; j < szMatrix; ++j) {
                if (j > i) {
                    matrixA[i * szMatrix + j] = matrixA[j * szMatrix + i] = i * szMatrix + j + 1;
                }
                matrixRotate[i * szMatrix + j] = matrixRotateT[i * szMatrix + j] =
                    matrixB[i * szMatrix + j] = matrixC[i * szMatrix + j] = 0;
            }
        }
        

        EPSk = Off(matrixA, cnt, kernelComputeNorm, command_queue, context, global_work_size);
        cout << "# " << EPSk << endl;
        while (EPSk >= EPS) {
            for (cl_uint k = 0; k < szMatrix; ++k) {
                matrixRotate[k * szMatrix + k] = 1;
            }
            cl_float* maxElem   = (cl_float*)malloc(sizeof(cl_float) * szMatrix);
            cl_int* maxPos      = (cl_int*)malloc(sizeof(cl_int) * max(1u * 2, szMatrix));
            
            cl_mem buffMaxPos       = NULL;
            cl_mem buffMaxElem      = NULL;

            buffMatrixA  = clCreateBuffer(context, CL_MEM_READ_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);    
            buffMaxElem  = clCreateBuffer(context, CL_MEM_READ_WRITE, szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);    
            buffMaxPos  = clCreateBuffer(context, CL_MEM_READ_WRITE, max(1u * 2, szMatrix) * sizeof(cl_int), NULL, &err);
            checkRet(err, __LINE__);    
            err = clEnqueueWriteBuffer(command_queue, buffMatrixA, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixA, 0, NULL, NULL);        
            checkRet(err, __LINE__);
            err = clSetKernelArg(kernelFindMaxAbs, 0, sizeof(cl_mem), (void *)&buffMatrixA);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelFindMaxAbs, 1, sizeof(cl_mem), (void *)&buffMaxElem);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelFindMaxAbs, 2, sizeof(cl_mem), (void *)&buffMaxPos);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelFindMaxAbs, 3, sizeof(cl_uint), (void *)&szMatrix);
            checkRet(err, __LINE__);                      
            err = clEnqueueNDRangeKernel(command_queue, kernelFindMaxAbs, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
            checkRet(err, __LINE__);
            err = clEnqueueReadBuffer(command_queue, buffMaxPos, CL_TRUE, 0, max(1u * 2, szMatrix) * sizeof(cl_int), maxPos, 0, NULL, NULL);
            checkRet(err, __LINE__);            
            clReleaseMemObject(buffMatrixA);    
            clReleaseMemObject(buffMaxElem);
            clReleaseMemObject(buffMaxPos);

            cl_int i = maxPos[0], j = maxPos[1];
            free(maxElem);
            free(maxPos);
            cout << "# " << i << " | " << j << endl;
            cl_float angel = fabs(matrixA[i * szMatrix + i] - matrixA[j * szMatrix + j])
                < 0.001 * EPS ? M_PI / 4 :
                (0.5 * atan(2.0 * matrixA[i * szMatrix + j]
                    / (matrixA[i * szMatrix + i] - matrixA[j * szMatrix + j])));
            matrixRotate[i * szMatrix + i] = matrixRotate[j * szMatrix + j] = cos(angel);
            matrixRotate[i * szMatrix + j] = -(matrixRotate[j * szMatrix + i] = sin(angel));

            for (size_t i = 0; i < szMatrix; ++i) {
                for (size_t j = 0; j < szMatrix; ++j) {
                    cout << matrixRotate[i * szMatrix + j] << " ";
                }
                cout << endl;
            }


            buffMatrixA  = clCreateBuffer(context, CL_MEM_READ_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);
            buffMatrixB  = clCreateBuffer(context, CL_MEM_WRITE_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);
            err = clEnqueueWriteBuffer(command_queue, buffMatrixA, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixRotate, 0, NULL, NULL);        
            checkRet(err, __LINE__);
            err = clSetKernelArg(kernelTranspose, 0, sizeof(cl_mem), (void *)&buffMatrixA);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelTranspose, 1, sizeof(cl_mem), (void *)&buffMatrixB);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelTranspose, 2, sizeof(cl_uint), (void *)&szMatrix);
            checkRet(err, __LINE__);                        
            err = clEnqueueNDRangeKernel(command_queue, kernelTranspose, 1, NULL, global_work_size, NULL, 0, NULL, NULL);
            checkRet(err, __LINE__);
            err = clEnqueueReadBuffer(command_queue, buffMatrixB, CL_TRUE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixRotateT, 0, NULL, NULL);
            checkRet(err, __LINE__);
            clReleaseMemObject(buffMatrixA);
            clReleaseMemObject(buffMatrixB);

            for (size_t i = 0; i < szMatrix; ++i) {
                for (size_t j = 0; j < szMatrix; ++j) {
                    cout << matrixRotateT[i * szMatrix + j] << " ";
                }
                cout << endl;
            }


            buffMatrixA  = clCreateBuffer(context, CL_MEM_READ_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);
            buffMatrixB  = clCreateBuffer(context, CL_MEM_READ_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);            
            buffMatrixC  = clCreateBuffer(context, CL_MEM_WRITE_ONLY, szMatrix * szMatrix * sizeof(cl_float), NULL, &err);
            checkRet(err, __LINE__);
            err = clEnqueueWriteBuffer(command_queue, buffMatrixA, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixRotateT, 0, NULL, NULL);        
            checkRet(err, __LINE__);
            err = clEnqueueWriteBuffer(command_queue, buffMatrixB, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixA, 0, NULL, NULL);        
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelMultMatrix, 0, sizeof(cl_mem), (void *)&buffMatrixA);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelMultMatrix, 1, sizeof(cl_mem), (void *)&buffMatrixB);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelMultMatrix, 2, sizeof(cl_mem), (void *)&buffMatrixC);
            checkRet(err, __LINE__);                        
            err = clSetKernelArg(kernelMultMatrix, 3, sizeof(cl_uint), (void *)&szMatrix);
            checkRet(err, __LINE__);                                    
            err = clEnqueueNDRangeKernel(command_queue, kernelMultMatrix, 2, NULL, global_work_size, NULL, 0, NULL, NULL);
            checkRet(err, __LINE__);
            err = clEnqueueReadBuffer(command_queue, buffMatrixC, CL_TRUE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixC, 0, NULL, NULL);
            checkRet(err, __LINE__);
            
            err = clEnqueueWriteBuffer(command_queue, buffMatrixA, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixC, 0, NULL, NULL);        
            checkRet(err, __LINE__);
            err = clEnqueueWriteBuffer(command_queue, buffMatrixB, CL_FALSE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixRotate, 0, NULL, NULL);        
            err = clSetKernelArg(kernelMultMatrix, 0, sizeof(cl_mem), (void *)&buffMatrixA);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelMultMatrix, 1, sizeof(cl_mem), (void *)&buffMatrixB);
            checkRet(err, __LINE__);            
            err = clSetKernelArg(kernelMultMatrix, 2, sizeof(cl_mem), (void *)&buffMatrixC);
            checkRet(err, __LINE__);                        
            err = clSetKernelArg(kernelMultMatrix, 3, sizeof(cl_uint), (void *)&szMatrix);
            checkRet(err, __LINE__);                                    
            err = clEnqueueNDRangeKernel(command_queue, kernelMultMatrix, 2, NULL, global_work_size, NULL, 0, NULL, NULL);
            checkRet(err, __LINE__);
            err = clEnqueueReadBuffer(command_queue, buffMatrixC, CL_TRUE, 0, 
                szMatrix * szMatrix * sizeof(cl_float), matrixA, 0, NULL, NULL);
            checkRet(err, __LINE__);

            clReleaseMemObject(buffMatrixA);
            clReleaseMemObject(buffMatrixB);
            clReleaseMemObject(buffMatrixC);

            EPSk = Off(matrixA, cnt, kernelComputeNorm, command_queue, context, global_work_size);
            cout << "#: " << EPSk << endl;
            for (size_t i = 0; i < szMatrix; ++i) {
                for (size_t j = 0; j < szMatrix; ++j) {
                    cout << matrixA[i * szMatrix + j] << " ";
                }
                cout << endl;
            }

        }


        free(matrixA);
        free(matrixB);
        free(matrixC);
        free(matrixRotate);
        free(matrixRotateT);
    }
    out.close();   
    free(source_str);
    free(devices);
    free(platforms);
    clReleaseContext(context);
    clReleaseKernel(kernelMultMatrix);
    clReleaseProgram(program);
    clReleaseCommandQueue(command_queue);     
}

void RotateMethodCPU() {
    string plotFileCPU = "plotDataCPU";
    ofstream out(plotFileCPU, ios::out);
    for (size_t cnt = 10; cnt < 800; cnt += 50) {
        //size_t cnt = atoi(argv[1]);
        size_t szMatrix = cnt;
        cl_float*   matrixA = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        cl_float*   matrixB = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        cl_float*   matrixC = (cl_float*)malloc(sizeof(cl_float) * szMatrix * szMatrix);
        for (size_t i = 0; i < szMatrix; ++i) {
            for (size_t j = i; j < szMatrix; ++j) {
                matrixA[i * szMatrix + j] = matrixA[j * szMatrix + i] = i * szMatrix + j + 1;
                matrixB[i * szMatrix + j] = matrixB[j * szMatrix + i] = i * szMatrix + j + 1;
            }
        }
        auto start = chrono::system_clock::now();
        for (size_t i = 0; i < szMatrix; ++i) {
            for (size_t j = 0; j < szMatrix; ++j) {
                cl_float temp = 0;
                for (size_t k = 0; k < szMatrix; ++k) {
                    temp += matrixA[i * szMatrix + k] * matrixB[k * szMatrix + j];    
                }
                matrixC[i * szMatrix + j] = temp;
            }
        }
        auto end = chrono::system_clock::now();
        auto elapsed_ns = chrono::duration_cast<chrono::milliseconds>(end - start);
        out << cnt << " " << elapsed_ns.count() << endl;
    }
    out.close();
}

int main(int argc, char* argv[]) {
    RotateMethodGPU();
    //RotateMethodCPU();
    return 0;
}
