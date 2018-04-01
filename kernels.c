/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/* 
 * Please fill in the following team struct 
 */
team_t team = {
"bcalabre + cemmel",              /* Team name */

    "Brent Calabresi",     /* First member full name */
    "bcalabre@u.rochester.edu",  /* First member email address */

    "Clay Emmel",                   /* Second member full name (leave blank if none) */
    "cemmel@u.rochester.edu"                    /* Second member email addr (leave blank if none) */
};
/***************
 * ROTATE KERNEL
 ***************/

/******************************************************
 * Your different versions of the rotate kernel go here
 ******************************************************/

/* 
 * naive_rotate - The naive baseline version of rotate 
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(dim-1-j, i, dim)] = src[RIDX(i, j, dim)];
}

/* 
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst) 
{
	int i,j,k ;
    	for (i = 0; i < dim; i = i + 16){
        	for (j = 0; j < dim; j++){
            		for(k = 0; k < 16; k++){
                   		dst[RIDX(dim-1-j, i+k, dim)] = src[RIDX(i+k, j, dim)];

             		}

        	}

    	}

}

/*********************************************************************
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_rotate_functions() 
{
    add_rotate_function(&naive_rotate, naive_rotate_descr);   
    add_rotate_function(&rotate, rotate_descr);   
    /* ... Register additional test functions here */
}


/***************
 * SMOOTH KERNEL
 **************/

/***************************************************************
 * Various typedefs and helper functions for the smooth function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
static int min(int a, int b) { return (a < b ? a : b); }
static int max(int a, int b) { return (a > b ? a : b); }

/* 
 * initialize_pixel_sum - Initializes all fields of sum to 0 
 */
static void initialize_pixel_sum(pixel_sum *sum) 
{
    sum->red = sum->green = sum->blue = 0;
    sum->num = 0;
    return;
}

/* 
 * accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum 
 */
static void accumulate_sum(pixel_sum *sum, pixel p) 
{
    sum->red += (int) p.red;
    sum->green += (int) p.green;
    sum->blue += (int) p.blue;
    sum->num++;
    return;
}

/* 
 * assign_sum_to_pixel - Computes averaged pixel value in current_pixel 
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum) 
{
    current_pixel->red = (unsigned short) (sum.red/sum.num);
    current_pixel->green = (unsigned short) (sum.green/sum.num);
    current_pixel->blue = (unsigned short) (sum.blue/sum.num);
    return;
}

/* 
 * avg - Returns averaged pixel value at (i,j) 
 */
static pixel avg(int dim, int i, int j, pixel *src) 
{
    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;

    initialize_pixel_sum(&sum);
    for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) 
	for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) 
	    accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);

    assign_sum_to_pixel(&current_pixel, sum);
    return current_pixel;
}

/******************************************************
 * Your different versions of the smooth kernel go here
 ******************************************************/

/*
 * naive_smooth - The naive baseline version of smooth 
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst) 
{
    int i, j;

    for (i = 0; i < dim; i++)
	for (j = 0; j < dim; j++)
	    dst[RIDX(i, j, dim)] = avg(dim, i, j, src);
}

char smooth_descr_Initial[] = "smooth: Initial Attempt";
void smooth_Initial(int dim, pixel *src, pixel *dst)
{
    int j = 0;
    int i = 0;

    //Loop unroll x16 - appears to be rather.. ineffective :(
    for (j = 0; j < dim; j+=16)
    {

        for (i = 0;  i< dim; i++)
        {

        dst[RIDX(j, i, dim)] = avg(dim, j, i, src);
        dst[RIDX(j+1, i, dim)] = avg(dim, j+1, i, src);
        dst[RIDX(j+2, i, dim)] = avg(dim, j+2, i, src);
        dst[RIDX(j+3, i, dim)] = avg(dim, j+3, i, src);
        dst[RIDX(j+4, i, dim)] = avg(dim, j+4, i, src);
        dst[RIDX(j+5, i, dim)] = avg(dim, j+5, i, src);
        dst[RIDX(j+6, i, dim)] = avg(dim, j+6, i, src);
        dst[RIDX(j+7, i, dim)] = avg(dim, j+7, i, src);
        dst[RIDX(j+8, i, dim)] = avg(dim, j+8, i, src);
        dst[RIDX(j+9, i, dim)] = avg(dim, j+9, i, src);
        dst[RIDX(j+10, i, dim)] = avg(dim, j+10, i, src);
        dst[RIDX(j+11, i, dim)] = avg(dim, j+11, i, src);
        dst[RIDX(j+12, i, dim)] = avg(dim, j+12, i, src);
        dst[RIDX(j+13, i, dim)] = avg(dim, j+13, i, src);
        dst[RIDX(j+14, i, dim)] = avg(dim, j+14, i, src);
        dst[RIDX(j+15, i, dim)] = avg(dim, j+15, i, src);
        }
    }
}


char smooth_descr_Best[] = "smooth: Best working version";
void smooth_Best(int dim, pixel *src, pixel *dst) 
{
int i, j, k, p;

    //Split the image into zones, top left corner, top right, bottom left, etc.

	
    /*CORNERS ----------------------------------------------------------------
    *
    * To get the avg of each corner we add the 2x2 matricies found in
    * each corner, then store the avg in the corner-most pixel
    */
   
    //Top Left Corner	
    dst[0].red =   (src[0].red + src[1].red + src[dim].red + src[dim+1].red) 	       >> 2;
    dst[0].blue =  (src[0].blue + src[1].blue + src[dim].blue + src[dim+1].blue)       >> 2;
    dst[0].green = (src[0].green + src[1].green + src[dim].green + src[dim + 1].green) >> 2;
    
    //Top Right Corner
    dst[dim-1].red =   (src[dim-1].red + src[dim-2].red + src[dim*2-1].red + src[dim*2-2].red) 	       >> 2;
    dst[dim-1].blue =  (src[dim-1].blue + src[dim-2].blue + src[dim*2-1].blue + src[dim*2-2].blue)     >> 2;
    dst[dim-1].green = (src[dim-1].green + src[dim-2].green + src[dim*2-1].green + src[dim*2-2].green) >> 2;
    
    //Bottom Left Corner
    dst[dim*(dim-1)].red =   (src[dim*(dim-1)].red + src[dim*(dim-1)+1].red + src[dim*(dim-2)].red + src[dim*(dim-2)+1].red)         >> 2;
    dst[dim*(dim-1)].blue =  (src[dim*(dim-1)].blue + src[dim*(dim-1)+1].blue + src[dim*(dim-2)].blue + src[dim*(dim-2)+1].blue)     >> 2;
    dst[dim*(dim-1)].green = (src[dim*(dim-1)].green + src[dim*(dim-1)+1].green + src[dim*(dim-2)].green + src[dim*(dim-2)+1].green) >> 2;
    
    //Bottom Right Corner
    dst[dim*dim-1].red = (src[dim*dim-1].red + src[dim*dim-2].red + src[dim*(dim-1)-1].red + src[dim*(dim-1)-2].red)           >> 2;
    dst[dim*dim-1].blue = (src[dim*dim-1].blue + src[dim*dim-2].blue + src[dim*(dim-1)-1].blue + src[dim*(dim-1)-2].blue)      >> 2;
    dst[dim*dim-1].green = (src[dim*dim-1].green + src[dim*dim-2].green + src[dim*(dim-1)-1].green + src[dim*(dim-1)-2].green) >> 2;
    

    /*EDGES ----------------------------------------------------------------
    *
    * Here, we iterate down the length of each edge and 
    * store the avg's for each point
    */

    //Top Side
    for (j = 1; j < dim - 1; j++) 
    {
        dst[j].red =   (src[j].red + src[j-1].red + src[j+1].red + src[j+dim].red + src[j+1+dim].red + src[j-1+dim].red) 		/ 6;
        dst[j].blue =  (src[j].blue + src[j-1].blue + src[j+1].blue + src[j+dim].blue + src[j+1+dim].blue + src[j-1+dim].blue) 		/ 6;
        dst[j].green = (src[j].green + src[j-1].green + src[j+1].green + src[j+dim].green + src[j+1+dim].green + src[j-1+dim].green) 	/ 6;
    }

    //Bottom Side
    for (j = dim * (dim - 1) + 1; j < dim * dim - 1; j++) 
    {
        dst[j].red =   (src[j].red + src[j-1].red + src[j+1].red + src[j-dim].red + src[j+1-dim].red + src[j-1-dim].red) 		/ 6;
        dst[j].blue =  (src[j].blue + src[j-1].blue + src[j+1].blue + src[j-dim].blue + src[j+1-dim].blue + src[j-1-dim].blue) 		/ 6;
        dst[j].green = (src[j].green + src[j-1].green + src[j+1].green + src[j-dim].green + src[j+1-dim].green + src[j-1-dim].green) 	/ 6;
    }
    //Left Hand Side
    for (j = dim; j < dim * (dim - 1); j += dim) 
    {
        dst[j].red =   (src[j].red + src[j-dim].red + src[j+1].red + src[j+dim].red + src[j+1+dim].red + src[j-dim+1].red) 		/ 6;
        dst[j].blue =  (src[j].blue + src[j-dim].blue + src[j+1].blue + src[j+dim].blue + src[j+1+dim].blue + src[j-dim+1].blue) 	/ 6;
        dst[j].green = (src[j].green + src[j-dim].green + src[j+1].green + src[j+dim].green + src[j+1+dim].green + src[j-dim+1].green) 	/ 6;
    }
    //Right Hand Side
    for (j = dim + dim - 1; j < dim * dim - 1; j += dim) 
    {
        dst[j].red =   (src[j].red + src[j-1].red + src[j-dim].red + src[j+dim].red + src[j-dim-1].red + src[j-1+dim].red) 		/ 6;
        dst[j].blue =  (src[j].blue + src[j-1].blue + src[j-dim].blue + src[j+dim].blue + src[j-dim-1].blue + src[j-1+dim].blue) 	/ 6;
        dst[j].green = (src[j].green + src[j-1].green + src[j-dim].green + src[j+dim].green + src[j-dim-1].green + src[j-1+dim].green) 	/ 6;
    }

   /*MIDDLE ----------------------------------
   * For the inner part of the matrix, we hard code the averages, very similar to the original naive system
   * But, since we handled the outer edges and corners already, less work is required and therefore less CPE's!
   */
    p = dim;
    for (i = 1; i < dim - 1; i++)
    {
        for (j = 1; j < dim - 1; j++)
        {
            k = p + j;
            dst[k].red =   (src[k].red + src[k-1].red + src[k+1].red + src[k-dim].red + src[k-dim-1].red + src[k-dim+1].red + src[k+dim].red + src[k+dim+1].red + src[k+dim-1].red) 			/ 9;
            dst[k].green = (src[k].green + src[k-1].green + src[k+1].green + src[k-dim].green + src[k-dim-1].green + src[k-dim+1].green + src[k+dim].green + src[k+dim+1].green + src[k+dim-1].green) 	/ 9;
            dst[k].blue =  (src[k].blue + src[k-1].blue + src[k+1].blue + src[k-dim].blue + src[k-dim-1].blue + src[k-dim+1].blue + src[k+dim].blue + src[k+dim+1].blue + src[k+dim-1].blue) 		/ 9;
        }
        p += dim;
    }
}


/*
 * smooth - Your current working version of smooth. 
 * IMPORTANT: This is the version you will be graded on
 */
char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst) 
{
 smooth_Best(dim, src, dst);
}

/********************************************************************* 
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_smooth_functions() {
    add_smooth_function(&smooth, smooth_descr);
    add_smooth_function(&naive_smooth, naive_smooth_descr);
    add_smooth_function(&smooth_Best,smooth_descr_Best);
    add_smooth_function(&smooth_Initial,smooth_descr_Initial);
    /* ... Register additional test functions here */
}

