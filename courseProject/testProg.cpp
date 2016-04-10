#include <CL/cl.h>
#include <iostream>
#include <cstdio>

#define MAX_SOURCE_SIZE (0x100000)

#define DEBUG 1

void checkRet(const cl_int err, const int raw) {
    if (err != CL_SUCCESS) {
        std::cerr << "Trouble: " << err << " | " << raw << std::endl;
        exit(1);
    }
}


int main() {
    cl_int              ret;
    cl_platform_id      platform_id;
    cl_uint             ret_num_platforms;
    cl_device_id        device_id;
    cl_uint             ret_num_devices;
    cl_context          context;
    cl_command_queue    command_queue;

    /* получить доступные платформы */
    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    checkRet(ret, __LINE__);
    /* получить доступные устройства */
    ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
    checkRet(ret, __LINE__);
    /* создать контекст */
    context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
    checkRet(ret, __LINE__);
    /* создаем команду */
    command_queue = clCreateCommandQueue(context, device_id, 0, &ret);    
    checkRet(ret, __LINE__);


    cl_program          program     = NULL;
    cl_kernel           kernel      = NULL;

    FILE* fp;
    const char* fileName = "./test.cl";
    char *source_str;
    size_t source_size;


    try {
        fp = fopen(fileName, "r");
        if (!fp) {
            std::cerr << "Failed to load kernel." << std::endl;
            exit(1);
        }
#ifdef DEBUG
        std::cout << "!" << std::endl;
#endif
        source_str = (char *)malloc(MAX_SOURCE_SIZE);
        source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
        fclose(fp);
    }
    catch (int a) {
        printf ("%d", a);
    }
    /* создать бинарник из кода программы */
    program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
    /* скомпилировать программу */
    checkRet(ret, __LINE__);
    ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    std::cout << ret << std::endl;
    checkRet(ret, __LINE__);

    /*  создать кернел */
    kernel = clCreateKernel(program, "test", &ret);
    checkRet(ret, __LINE__);
    
    size_t global_work_size = 10;
    cl_mem      memobj      = NULL;
    cl_int      memLength   = 10;
    cl_double   **mem       = (cl_double **)malloc(sizeof(cl_double *) * memLength);
    for (int i = 0; i < memLength; ++i) {
        mem[i] = (cl_double *)malloc(sizeof(cl_double) * memLength);
        for (int j = 0; j < memLength; ++j)
            mem[i][j] = j;
    }

    for (int i = 0; i < memLength; ++i) {
        /* создать буфер */
        memobj  = clCreateBuffer(context, CL_MEM_READ_WRITE, global_work_size * sizeof(cl_float), NULL, &ret);

        /* записать данные в буфер */
        ret = clEnqueueWriteBuffer(command_queue, memobj, CL_TRUE, 0, global_work_size * sizeof(cl_float), mem[i], 0, NULL, NULL);        

        if (!!ret) {
            std::cerr << "Beda0" << std::endl;
            exit(1);
        }
        /* устанавливаем параметр (2-ой параметр отвечает за порядковый номер 4-ого аргумента в ф-ии __kernel */
        ret = clSetKernelArg(kernel, 0, sizeof(cl_double*), (void *)&mem[i]);
        if (!!ret) {
            std::cerr << "Beda1" << std::endl;
            exit(2);
        }
//        ret = clSetKernelArg(kernel, 1, sizeof(cl_int), (void *)&memLength);
//        if (!!ret) {
//            std::cerr << "Beda2" << std::endl;
//        }

#ifdef DEBUG
        std::cout << "#" << std::endl;
#endif

        /* выполнить кернел */
        ret = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);

#ifdef DEBUG
        std::cout << "&" << std::endl;
#endif
        /* считать данные из буфера */
        ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0, memLength * sizeof(cl_float), mem[i], 0, NULL, NULL);
        for (int j = 0; j < memLength; ++j)
            std::cout << mem[i][j] << " ";
        std::cout << std::endl;
    }

    for (int i = 0; i < memLength; ++i) {
        for (int j = 0; j < memLength; ++j)
            std::cout << mem[i][j] << " ";
        std::cout << std::endl;
    }
    return 0;
}
