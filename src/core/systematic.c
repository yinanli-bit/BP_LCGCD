#include "bplcgcd/gf2.h"

#include <stdlib.h>
#include <string.h>

void binary_systematic_free(BinarySystematic *form)
{
    if (form) {
        binary_matrix_free(&form->reduced);
        binary_matrix_free(&form->row_transform);
        free(form->sys_to_orig);
        free(form->orig_to_sys);
        *form = (BinarySystematic){0};
    }
}

static void swap_rows(BinaryMatrix *matrix, size_t a, size_t b)
{
    for (size_t j = 0; j < matrix->cols; ++j) {
        uint8_t value = matrix->data[a * matrix->cols + j];
        matrix->data[a * matrix->cols + j] = matrix->data[b * matrix->cols + j];
        matrix->data[b * matrix->cols + j] = value;
    }
}

//得到秩、D、U 和列映射
Gf2Status binary_matrix_systematize(const BinaryMatrix *matrix,BinarySystematic *form)
{
    if (!matrix || !matrix->data || !matrix->rows || !matrix->cols ||
        matrix->rows > SIZE_MAX / matrix->cols ||
        matrix->cols > SIZE_MAX / sizeof(size_t) || !form || form->rank ||
        form->reduced.data || form->reduced.rows || form->reduced.cols ||
        form->row_transform.data || form->row_transform.rows ||
        form->row_transform.cols || form->sys_to_orig || form->orig_to_sys) {
        return GF2_INVALID_ARGUMENT;
    }
    size_t rows = matrix->rows, cols = matrix->cols;
    for (size_t i = 0; i < rows * cols; ++i) {
        if (matrix->data[i] > 1) return GF2_INVALID_ARGUMENT;
    }
    BinarySystematic built = {0}; //结构体全部初始化为0
    Gf2Status status = binary_matrix_init(&built.reduced, rows, cols);
    if (status == GF2_OK) {
        status = binary_matrix_init(&built.row_transform, rows, rows);
    }
    if (status != GF2_OK) {
        binary_systematic_free(&built);
        return status;
    }
    built.sys_to_orig = malloc(cols * sizeof(size_t));
    built.orig_to_sys = malloc(cols * sizeof(size_t));
    if (!built.sys_to_orig || !built.orig_to_sys) {
        binary_systematic_free(&built);
        return GF2_ALLOCATION_FAILED;
    }
    uint8_t *a = built.reduced.data, *u = built.row_transform.data;
    memcpy(a, matrix->data, rows * cols);
    for (size_t i = 0; i < rows; ++i) u[i * rows + i] = 1;
    for (size_t j = 0; j < cols; ++j) built.sys_to_orig[j] = j;
    for (; built.rank < rows && built.rank < cols; ++built.rank) {
        //第一次寻找主元时build.rank=0
        size_t k = built.rank, pivot_col, pivot_row = rows;
        /* Deterministic: first available column, then first available row. */
        for (pivot_col = k; pivot_col < cols; ++pivot_col) {
            for (pivot_row = k; pivot_row < rows; ++pivot_row) {
                if (a[pivot_row * cols + pivot_col]) break;
            }
            if (pivot_row < rows) break;
        }
         //找到主元坐标(pivot_row,pivot_col)
        
        printf("=====第%d次寻找=====\n",k);
        printf("(%d,%d)\n",pivot_row,pivot_col);
        
        if (pivot_col == cols) break;
        swap_rows(&built.reduced, k, pivot_row);
        swap_rows(&built.row_transform, k, pivot_row);
        for (size_t i = 0; i < rows; ++i) {
            uint8_t value = a[i * cols + k];
            a[i * cols + k] = a[i * cols + pivot_col];
            a[i * cols + pivot_col] = value;
        }//交换列
        //记录列映射关系，原始-->系统
        size_t original = built.sys_to_orig[k]; 
        built.sys_to_orig[k] = built.sys_to_orig[pivot_col];
        built.sys_to_orig[pivot_col] = original;
        //printf("%s\n",built.sys_to_orig);
        
        printf("sys_to_orig:\n");
        for(int i=0;i<cols;++i){
            printf("%d ",built.sys_to_orig[i]);
        }
        printf("\n");
        

        for (size_t i = 0; i < rows; ++i) {
            if (i == k || !a[i * cols + k]) continue;
            for (size_t j = 0; j < cols; ++j) a[i * cols + j] ^= a[k * cols + j];
            for (size_t j = 0; j < rows; ++j) u[i * rows + j] ^= u[k * rows + j];//记录行消元模式
        }//消元
        
        printf("消元矩阵a:\n");
        for(int i=0;i<rows;++i){
            for(int j=0;j<cols;++j){
                printf("%d ",a[i*cols+j]);
            }
            printf("\n");
        }
        printf("消元矩阵u:\n");
        for(int i=0;i<rows;++i){
            for(int j=0;j<rows;++j){
                printf("%d ",u[i*rows+j]);
            }
            printf("\n");
        }
        
        
    } 
    for (size_t j = 0; j < cols; ++j) {
        built.orig_to_sys[built.sys_to_orig[j]] = j;
    }//将系统顺序恢复成原始顺序，orig_to_sys[0]=2：原始第0列在系统化矩阵中的第2列，错误向量也要跟着变
    //printf("%s\n",built.orig_to_sys);
    
    printf("orig_to_sys:\n");
    for(int i=0;i<cols;++i){    
        printf("%d ",built.orig_to_sys[i]);
    }
    printf("\n");
    
    *form = built;
    return GF2_OK;
}

static int ready(const BinarySystematic *form)
{
    return form && form->reduced.data && form->row_transform.data &&
           form->sys_to_orig && form->orig_to_sys &&
           form->rank <= form->reduced.rows && form->rank <= form->reduced.cols;
}

static Gf2Status permute(const BinarySystematic *form, const BitVector *input,
                         BitVector *result, int inverse)//reverse=1表示从系统化矩阵映射回原始矩阵，inverse=0表示从原始矩阵映射到系统化矩阵
{
    if (!ready(form) || !input || !result || !input->data || !result->data ||
        result->data == input->data || result->data == form->reduced.data ||
        result->data == form->row_transform.data) {
        return GF2_INVALID_ARGUMENT;
    }
    size_t cols = form->reduced.cols;
    if (input->length != cols || result->length != cols) {
        return GF2_DIMENSION_MISMATCH;
    }
    for (size_t j = 0; j < cols; ++j) {
        if (input->data[j] > 1) return GF2_INVALID_ARGUMENT;
    }
    const size_t *map = inverse ? form->orig_to_sys : form->sys_to_orig;
    for (size_t j = 0; j < cols; ++j) result->data[j] = input->data[map[j]];
    return GF2_OK;
}

Gf2Status binary_systematic_to_system(const BinarySystematic *form,
                                    const BitVector *input, BitVector *result)
{
    return permute(form, input, result, 0);
}

Gf2Status binary_systematic_to_original(const BinarySystematic *form,
                                      const BitVector *input, BitVector *result)
{
    return permute(form, input, result, 1);
}

//syndrome变换：s_sys = U * s_orig；同时检查一致性：s_sys的后面部分必须为0
Gf2Status binary_systematic_transform_syndrome(const BinarySystematic *form,
                                            const BitVector *input,
                                            BitVector *result)
{
    if (!ready(form) || !result || result->data == form->reduced.data) {
        return GF2_INVALID_ARGUMENT;
    }
    Gf2Status status = binary_matrix_multiply_vector(&form->row_transform,input, result);
    if (status != GF2_OK) return status;
    for (size_t i = form->rank; i < result->length; ++i) {
        if (result->data[i]) return GF2_INCONSISTENT_SYNDROME;
    }
    
    printf("s_sys:\n");
    for(int i=0;i<result->length;++i){
        printf("%d ",result->data[i]);
    }
    printf("\n");
    

    return GF2_OK;
}
