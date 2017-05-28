/**
 * @file    radix_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 *
 * Copyright (c) 2016 Georgia Institute of Technology. All Rights Reserved.
 */

/*
 * TODO: implement your radix sort solution in this file
 */

#include <mpi.h>

// returns the value of the digit starting at offset `offset` and containing `k` bits
#define GET_DIGIT(key, k, offset) ((key) >> (offset)) & ((1 << (k)) - 1)


/**
 * @brief   Parallel distributed radix sort.
 *
 * This function sorts the distributed input range [begin, end)
 * via lowest-significant-byte-first radix sort.
 *
 * This function will sort elements of type `T`, by the key of type `unsigned int`
 * which is returned by the key-access function `key_func`.
 *
 * The MPI datatype for the templated (general) type `T` has to be passed
 * via the `dt` parameter.
 *
 * @param begin         A pointer to the first element in the local range to be sorted.
 * @param end           A pointer to the end of the range to be sorted. This
 *                      pointer points one past the last element, such that the
 *                      total number of elements is given by `end - begin`.
 * @param key_func      A function with signature: `unsigned int (const T&)`.
 *                      This function returns the key of each element, which is
 *                      used for sorting.
 * @param dt            The MPI_Datatype which represents the type `T`. This
 *                      is used whenever elements of type `T` are communicated
 *                      via MPI.
 * @param comm          The communicator on which the sorting happens.
 *                      NOTE: this is not necessarily MPI_COMM_WORLD. Call
 *                            all MPI functions with this communicator and
 *                            NOT with MPI_COMM_WORLD.
 */
template <typename T>
void radix_sort(T* begin, T* end, unsigned int (*key_func)(const T&), MPI_Datatype dt, MPI_Comm comm, unsigned int k = 16) {
    // get comm rank and size
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // The number of elements per processor: n/p
    size_t np = end - begin;

    // the number of histogram buckets = 2^k
    unsigned int num_buckets = 1 << k;

    for (unsigned int d = 0; d < 8*sizeof(unsigned int); d += k) {
        // TODO:
        // 1.) create histogram and sort via bucketing (~ counting sort)
        std::vector<unsigned int> histogram(num_buckets);
        std::vector<T> crops(np);
        //gets histogram
        for (unsigned int a = 0; a < np; a++) {
            T address = *(begin + a);
            unsigned int value = key_func(address);
            crops[a] = address;
            unsigned int digit = GET_DIGIT(value, k, d);
            histogram[digit] += 1;
        }
        std::vector<unsigned int> positions(num_buckets);
        positions[0] = 0;
        //calculate starting positions for each digit in the local order
        for (unsigned int b = 1; b < num_buckets; b++) {
             unsigned int firstVal = positions[b - 1];
             unsigned int secVal = histogram[b - 1];
             positions[b] = firstVal + secVal;
        }
        //calculate sorted local positions
        std::vector<unsigned int> sortedPos(np);
        for (unsigned int a = 0; a < np; a++) {
            T address = *(begin + a);
            unsigned int value = key_func(address);
            unsigned int digit = GET_DIGIT(value, k, d);
            sortedPos[a] = positions[digit];
            positions[digit] += 1;
        }
        //put elements in proper sorted order
        for (unsigned int a = 0; a < np; a++) {
            T* newAddress = begin + sortedPos[a];
            *newAddress = crops[a];
        }
        // 2.) get global histograms (P, G) via MPI_Exscan/MPI_Allreduce,...
        //std::vector<unsigned int> localG(np);
        //for (unsigned int a = 0; a < np; a++) {
        //      T address = *(begin + a);
        //      unsigned int value = key_func(address);
        //      unsigned int digit = GET_DIGIT(value, k, d*k);
        //      int totalSum = 0;
        //      if (digit == 0) {
        //        localG[a] = 0;
        //      } else {
        //	      for (unsigned int c = 0; c < digit; c++) {
        //         totalSum += histogram[c]; 
        //        }
        //      localG[a] = totalSum;
        //      }
        //construct local P histogram
        //std::vector<unsigned int> localP(np);
        //for (unsigned int a = 0; a < np; a++) {
        //      T address = *(begin + a);
        //      unsigned int value = key_func(address);
        //      unsigned int digit = GET_DIGIT(value, k, d*k);
        //      localP[a] = histogram[digit];
        //}
        //get local L values
        std::vector<unsigned int> localL(np);
        int currSum = 0;
        int currDigit = 0;
        for (unsigned int b = 0; b < np; b++) {
             T address = *(begin + b);
             unsigned int value = key_func(address);
             unsigned int digit = GET_DIGIT(value, k, d);
             if (digit != currDigit) {
                currSum = 0;
                currDigit = digit;
             }
             localL[b] = currSum;
             currSum += 1;
        }
        //use MPI_ALLreduce to get global G
        std::vector<unsigned int> globalG(num_buckets);
        MPI_Allreduce(&histogram[0], &globalG[0], num_buckets, MPI_UNSIGNED, MPI_SUM, comm);
        //use MPI_Exscan to get global P
        std::vector<unsigned int> globalP(num_buckets);
        MPI_Exscan(&histogram[0], &globalP[0], num_buckets, MPI_UNSIGNED, MPI_SUM, comm);
        std::vector<unsigned int> signedCounts(np);
        for (unsigned int z = 0; z < np; z++) {
            T address = *(begin + z);
            unsigned int value = key_func(address);
            unsigned int digit = GET_DIGIT(value, k, d);
            unsigned int theG = 0;
            for (int c = 0; c < digit; c++) {
                //printf("the current value is %d and the current value of digit is %d, totes is %d",value, digit, globalG[c]);
                theG += globalG[c];
            }
            unsigned int theP = globalP[digit];
            signedCounts[z] = localL[z] + theG + theP;
        }
        //calculate send_counts for all of the processors, how many elems to send
        std::vector<int> send_counts(p);
        std::vector<T> elements(np);
        for (unsigned int f = 0; f < np; f++) {
             int robo = (signedCounts[f]) / (np);
             send_counts[robo] += 1;
        }
        //calculate displacements
        std::vector<int> displacements(p);
        displacements[0] = 0;
        for (unsigned int e = 1; e < p; e++) {
            displacements[e] = displacements[e - 1] + send_counts[e - 1];
        }
        std::vector<int> recv_counts(p);
        //get recv_c counts
        MPI_Alltoall(&send_counts[0], 1, MPI_INT, &recv_counts[0], 1, MPI_INT, comm);
        std::vector<int> recvdisplacements(p);
        recvdisplacements[0] = 0;
        for (unsigned int f = 1; f < p; f++) {
            recvdisplacements[f] = recvdisplacements[f - 1] + recv_counts[f - 1];
        }
        // put elements in right global order
        MPI_Alltoallv(begin, &send_counts[0], &displacements[0], dt, &elements[0], &recv_counts[0], &recvdisplacements[0], dt, comm);
        //local sorting via bucketing (~ counting sort), same procedure as before
        //uses elements and sorts elems back into old begin
        std::vector<unsigned int> new_histogram(num_buckets);
        T* newBegin = &elements.at(0);
        for (unsigned int a = 0; a < np; a++) {
            T theStruct = *(newBegin + a);
            unsigned int value = key_func(theStruct);
            unsigned int digit = GET_DIGIT(value, k, d);
            new_histogram[digit] += 1;
        }
        std::vector<unsigned int> newPositions(num_buckets);
        newPositions[0] = 0;
        for (unsigned int b = 1; b < num_buckets; b++) {
             unsigned int firstVal = (unsigned int) newPositions[b - 1];
             unsigned int secVal = (unsigned int)new_histogram[b - 1];
             newPositions[b] = firstVal + secVal;
        }
        std::vector<unsigned int> sortedP(np);
        for (unsigned int a = 0; a < np; a++) {
            T address = *(newBegin + a);
            unsigned int value = key_func(address);
            unsigned int digit = GET_DIGIT(value, k, d);
            sortedP[a] = newPositions[digit];
            newPositions[digit] += 1;
        }     
        for (unsigned int a = 0; a < np; a++) {
            T* newAddress = begin + sortedP[a];          
            *newAddress = elements[a];
        }
    }
}

  
