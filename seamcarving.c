#include "seamcarving.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/*

struct rgb_img{
    uint8_t *raster;
    size_t height;
    size_t width;
};

*/

void calc_energy(struct rgb_img *im, struct rgb_img **grad)
{
    *grad = (struct rgb_img *)malloc(sizeof(struct rgb_img));
    (*grad)->height = im->height;
    (*grad)->width = im->width;
    (*grad)->raster = (uint8_t *)malloc(3 * (*grad)->height * (*grad)->width);
    for (int i = 0; i < im->height; i++)
    {
        for (int j = 0; j < im->width; j++)
        {
            int rx = (int) (get_pixel(im, i, ((j-1) + im->width) % im->width, 0) - get_pixel(im, i, (j+1) % im->width, 0));
            rx = rx*rx;
            int gx = (int) (get_pixel(im, i, ((j-1) + im->width) % im->width, 1) - get_pixel(im, i, (j+1) % im->width, 1));
            gx = gx*gx;
            int bx = (int) (get_pixel(im, i, ((j-1) + im->width) % im->width, 2) - get_pixel(im, i, (j+1) % im->width, 2));
            bx = bx*bx; 
            int ry = (int) (get_pixel(im, ((i-1) + im->height) % im->height, j, 0) - get_pixel(im, (i+1) % im->height, j, 0));
            ry = ry*ry; 
            int gy = (int) (get_pixel(im, ((i-1) + im->height) % im->height, j, 1) - get_pixel(im, (i+1) % im->height, j, 1));
            gy = gy*gy; 
            int by = (int) (get_pixel(im, ((i-1) + im->height) % im->height, j, 2) - get_pixel(im, (i+1) % im->height, j, 2));
            by = by*by;
            int delta = rx+gx+bx+ry+gy+by;
            delta = (int) (sqrt(delta)/10);
            uint8_t delta8 = (uint8_t) delta;
            set_pixel(*grad, i, j, delta8, delta8, delta8);
        }
    }
}
void dynamic_seam(struct rgb_img *grad, double **best_arr)
{
    *best_arr = (double *)malloc(sizeof(double) * grad->width * grad->height);
    for (int i = 0; i < grad->width; i++)
    {
        (*best_arr)[i] = (double) get_pixel(grad, 0, i, 0);
    }
    for (int j = 1; j < grad->height; j++)
    {
        double minimum = ((*best_arr)[(j-1)*grad->width] < (*best_arr)[(j-1)*grad->width + 1]) ? (*best_arr)[(j-1)*grad->width] : (*best_arr)[(j-1)*grad->width + 1];
        (*best_arr)[j * grad->width] = minimum + (double) get_pixel(grad, j, 0, 0);
        for (int k = 1; k < grad->width - 1; k++)
        {
            double minimum2 = ((*best_arr)[(j-1)*grad->width + k] < (*best_arr)[(j-1)*grad->width + 1 + k]) ? (*best_arr)[(j-1)*grad->width + k] : (*best_arr)[(j-1)*grad->width + 1 + k];
            minimum2 = ((*best_arr)[(j-1)*grad->width - 1 + k] < minimum2) ? (*best_arr)[(j-1)*grad->width - 1 + k] : minimum2;
            (*best_arr)[j * grad->width + k] = minimum2 + (double) get_pixel(grad, j, k, 0);
        }
        double minimum3 = ((*best_arr)[j*grad->width - 2] < (*best_arr)[j*grad->width - 1]) ? (*best_arr)[j*grad->width - 2] : (*best_arr)[j*grad->width - 1];
        (*best_arr)[(j+1) * grad->width - 1] = minimum3 + (double) get_pixel(grad, j, grad->width - 1, 0);
    }
    free(best_arr);
}
void recover_path(double *best, int height, int width, int **path)
{
    *path = (int *)malloc(sizeof(int) * height);
    
    double min_energy = best[(height - 1) * width];
    int min_index = 0;
    for (int i = 1; i < width; i++) 
    {
        if (best[(height - 1) * width + i] < min_energy) 
        {
            min_energy = best[(height - 1) * width + i];
            min_index = i;
        }
    }

    (*path)[height - 1] = min_index;
    for (int i = height - 2; i >= 0; i--) {
        double cur_energy = best[i * width + min_index];

        double left_energy = 63; // 63 is the maximum energy
        if (min_index > 0)
        {
            left_energy = best[i * width + min_index - 1];
        }

        double middle_energy = cur_energy;

        double right_energy = 63;
        if (min_index < width - 1)
        {
            right_energy = best[i * width + min_index + 1];
        }

        if (left_energy <= middle_energy && left_energy <= right_energy)
            min_index -= 1;
        else if (right_energy <= middle_energy && right_energy <= left_energy)
            min_index += 1;
        
        (*path)[i] = min_index;
    }
}
void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path)
{
    *dest = (struct rgb_img *)malloc(sizeof(struct rgb_img));
    (*dest)->height = src->height;
    (*dest)->width = src->width - 1;
    (*dest)->raster = (uint8_t *)malloc(3 * (*dest)->height * (*dest)->width);
    
    for (int i = 0; i < src->height; i++) {
        int dest_index = 0;
        for (int j = 0; j < src->width; j++) {
            if (j != path[i]) {
                (*dest)->raster[3 * (i * (*dest)->width + dest_index)] = src->raster[3 * (i * src->width + j)];
                (*dest)->raster[3 * (i * (*dest)->width + dest_index) + 1] = src->raster[3 * (i * src->width + j) + 1];
                (*dest)->raster[3 * (i * (*dest)->width + dest_index) + 2] = src->raster[3 * (i * src->width + j) + 2];
                dest_index++;
            }
        }
    }
}


int main()
{
    
    struct rgb_img *grad;
    struct rgb_img *im;
    read_in_img(&im, "6x5.bin");
    calc_energy(im, &grad);
    print_grad(grad);
    double *best_arr;
    dynamic_seam(grad, &best_arr);
    printf("%f\n", *(best_arr+17));
    int *path;
    recover_path(best_arr, 5, 6, &path);
    printf("%d\n", *(path+4));
    struct rgb_img *dest;
    remove_seam(im, &dest, path);
    write_img(dest, "new6x5.bin");

    struct rgb_img *im;
    struct rgb_img *cur_im;
    struct rgb_img *grad;
    double *best;
    int *path;

    read_in_img(&im, "HJoceanSmall.bin");
    
    for(int i = 0; i < 5; i++){
        printf("i = %d\n", i);
        calc_energy(im,  &grad);
        dynamic_seam(grad, &best);
        recover_path(best, grad->height, grad->width, &path);
        remove_seam(im, &cur_im, path);

        char filename[200];
        sprintf(filename, "img%d.bin", i);
        write_img(cur_im, filename);


        destroy_image(im);
        destroy_image(grad);
        free(best);
        free(path);
        im = cur_im;
    }
    destroy_image(im);
}
