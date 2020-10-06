//
//  Sort.c
//  OLLP
//
//  Created by 高文智 on 2020/8/19.
//  Copyright © 2020 高文智. All rights reserved.
//

#include "Sort.h"
#include <stdio.h>

void swap(double *a, OLLP_int i, OLLP_int j)
{
    OLLP_int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
}

void swap2(OLLP_int *a, OLLP_int i, OLLP_int j)
{
    OLLP_int tmp = a[i];
    a[i] = a[j];
    a[j] = tmp;
}

OLLP_int Partition(OLLP_int* idx, double* L, OLLP_int low, OLLP_int high)
{
    OLLP_int tmp = low;
    OLLP_int pt = L[low];

    while (low < high) {
        
        while (low < high && L[high] >= pt) -- high;
        while (low < high && L[low] <= pt) ++low;
        if (low < high) {
            swap(L, low, high);
            swap2(idx, low, high);
        }
    }
    
    swap(L, low, tmp);
    swap2(idx, low, tmp);
    
    return low;
}

void QSort(OLLP_int* idx, double *L, OLLP_int low, OLLP_int high)
{
    if (low < high) {
        OLLP_int pl = Partition(idx, L, low, high);
        QSort(idx,L,low,pl - 1);
        QSort(idx,L,pl + 1,high);
    }
}


