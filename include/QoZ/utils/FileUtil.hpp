//
// Created by Kai Zhao on 12/9/19.
//

#ifndef _SZ_FILE_UTIL
#define _SZ_FILE_UTIL

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <random>
#include <sstream>
#include <cstdio>
int RW_SCES=0;
int RW_TERR=1;
int RW_FERR=2;
typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;
namespace QoZ {

    template<typename Type>
    void readfile(const char *file, const size_t num, Type *data) {
        std::ifstream fin(file, std::ios::binary);
        if (!fin) {
            std::cout << " Error, Couldn't find the file: " << file << "\n";
            exit(0);
        }
        fin.seekg(0, std::ios::end);
        const size_t num_elements = fin.tellg() / sizeof(Type);
        assert(num_elements == num && "File size is not equals to the input setting");
        fin.seekg(0, std::ios::beg);
        fin.read(reinterpret_cast<char *>(data), num_elements * sizeof(Type));
        fin.close();
    }

    template<typename Type>
    std::unique_ptr<Type[]> readfile(const char *file, size_t &num) {
        std::ifstream fin(file, std::ios::binary);
        if (!fin) {
            std::cout << " Error, Couldn't find the file: " << file << std::endl;
            exit(0);
        }
        fin.seekg(0, std::ios::end);
        const size_t num_elements = fin.tellg() / sizeof(Type);
        fin.seekg(0, std::ios::beg);
//        auto data = QoZ::compat::make_unique<Type[]>(num_elements);
        auto data = std::make_unique<Type[]>(num_elements);
        fin.read(reinterpret_cast<char *>(&data[0]), num_elements * sizeof(Type));
        fin.close();
        num = num_elements;
        return data;
    }

    template<typename Type>
    void writefile(const char *file, Type *data, size_t num_elements) {
        std::ofstream fout(file, std::ios::binary);
        fout.write(reinterpret_cast<const char *>(&data[0]), num_elements * sizeof(Type));
        fout.close();
    }

    template<typename Type>
    void writeTextFile(const char *file, Type *data, size_t num_elements) {
        std::ofstream fout(file);
        if (fout.is_open()) {
            for (size_t i = 0; i < num_elements; i++) {
                fout << data[i] << std::endl;
            }
            fout.close();
        } else {
            std::cout << "Error, unable to open file for output: " << file << std::endl;
            exit(0);
        }
    }


    //temp for test
    void writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status)
    {
        FILE *pFile = fopen(tgtFilePath, "wb");
        if (pFile == NULL)
        {
            printf("Failed to open input file. 3\n");
            *status = RW_FERR;
            return;
        }
        
        fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
        fclose(pFile);
        *status = RW_SCES;
    }

    void writeDoubleData(double *data, size_t nbEle, char *tgtFilePath, int *status)
    {
        size_t i = 0;
        char s[64];
        FILE *pFile = fopen(tgtFilePath, "wb");
        if (pFile == NULL)
        {
            printf("Failed to open input file. 3\n");
            *status = RW_FERR;
            return;
        }
        
        for(i = 0;i<nbEle;i++)
        {
            sprintf(s,"%.20G\n",data[i]);
            fputs(s, pFile);
        }
        
        fclose(pFile);
        *status = RW_SCES;
    }

    void writeFloatData(float *data, size_t nbEle, char *tgtFilePath, int *status)
    {
        size_t i = 0;
        char s[64];
        FILE *pFile = fopen(tgtFilePath, "wb");
        if (pFile == NULL)
        {
            printf("Failed to open input file. 3\n");
            *status = RW_FERR;
            return;
        }
       
        for(i = 0;i<nbEle;i++)
        {
            //printf("i=%d\n",i);
            //printf("data[i]=%f\n",data[i]);
            sprintf(s,"%.30G\n",data[i]);
            fputs(s, pFile);
        }
        
        fclose(pFile);
        *status = RW_SCES;
    }



    void writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status)
    {
        size_t i = 0; 
        int state = RW_SCES;
        lfloat buf;
        unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
        for(i=0;i<nbEle;i++)
        {
            buf.value = data[i];
            bytes[i*4+0] = buf.byte[0];
            bytes[i*4+1] = buf.byte[1];
            bytes[i*4+2] = buf.byte[2];
            bytes[i*4+3] = buf.byte[3];                 
        }

        size_t byteLength = nbEle*sizeof(float);
        writeByteData(bytes, byteLength, tgtFilePath, &state);
        free(bytes);
        *status = state;
    }

    void writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status)
    {
        size_t i = 0, index = 0; 
        int state = RW_SCES;
        ldouble buf;
        unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
        for(i=0;i<nbEle;i++)
        {
            index = i*8;
            buf.value = data[i];
            bytes[index+0] = buf.byte[0];
            bytes[index+1] = buf.byte[1];
            bytes[index+2] = buf.byte[2];
            bytes[index+3] = buf.byte[3];
            bytes[index+4] = buf.byte[4];
            bytes[index+5] = buf.byte[5];
            bytes[index+6] = buf.byte[6];
            bytes[index+7] = buf.byte[7];
        }

        size_t byteLength = nbEle*sizeof(double);
        writeByteData(bytes, byteLength, tgtFilePath, &state);
        free(bytes);
        *status = state;
    }

}

#endif