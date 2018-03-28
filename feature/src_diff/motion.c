/**
 *
 *  Copyright 2016-2017 Netflix, Inc.
 *
 *     Licensed under the Apache License, Version 2.0 (the "License");
 *     you may not use this file except in compliance with the License.
 *     You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 *     Unless required by applicable law or agreed to in writing, software
 *     distributed under the License is distributed on an "AS IS" BASIS,
 *     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *     See the License for the specific language governing permissions and
 *     limitations under the License.
 *
 */

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "motion_options.h"
#include "common_christos/tools.h"
#include "common_christos/alloc.h"
#include "common_christos/file_io.h"
#include "common_christos/convolution.h"
#include "common_christos/convolution_internal.h"
#include "motion_tools.h"

#ifdef MOTION_OPT_SINGLE_PRECISION
    typedef float number_t;
    #define read_image_b  read_image_b2s
    #define read_image_w  read_image_w2s
    #define convolution_f32_c convolution_f32_c_s
    #define FILTER_5           FILTER_5_s

#else
    typedef double number_t;
    #define read_image_b  read_image_b2d
    #define read_image_w  read_image_w2d
    #define convolution_f32_c convolution_f32_c_d
    #define FILTER_5           FILTER_5_d
#endif

/**
 * Note: img1_stride and img2_stride are in terms of (sizeof(number_t) bytes)
 */

number_t frame_differ(const number_t *img_after, const number_t *img_before, number_t *img_diff, int width, int height, int img_stride)
{

for (int i = 0; i < height; ++i) {
	for (int j = 0; j < width; ++j) {
		number_t img_px_after = img_after[i * img_stride + j];
		number_t img_px_before = img_before[i * img_stride + j];
		number_t lol = img_after[i * img_stride + j] - img_before[i * img_stride + j];
		img_diff[i * img_stride + j] = lol;
		number_t img_px_diff = img_diff[i * img_stride + j];
		//printf("i: %f, i+1: %f, diff: %f\n", img_px_after, img_px_before, img_px_diff);
	}
}

}

number_t vmaf_image_sad_c(const number_t *img1, const number_t *img2, int width, int height, int img1_stride, int img2_stride, int frm_idx, int blk_x, int blk_y, int new_rows, int new_cols, int p_o, int Nlines_x, int Nlines_y, int left_over_x, int left_over_y)
{
    number_t accum = (number_t)0.0;
    number_t motion_blocks_row[Nlines_x * Nlines_y];
    int row_accum, col_accum;

    /*int left_over_y = width - (new_cols) * blk_y;
    int left_over_x = height - (new_rows) * blk_x;*/

    for (int i = 0; i < Nlines_x * Nlines_y; i++) { 
		motion_blocks_row[i] = 0; 
    }

    for (int block_start_x_case = 1; block_start_x_case < Nlines_x + 1; block_start_x_case++) {
        for (int block_start_y_case = 1; block_start_y_case < Nlines_y + 1; block_start_y_case++) {
	    int top_now = (block_start_x_case - 1) * (blk_x / p_o) + (left_over_x / 2);
	    int left_now = (block_start_y_case - 1) * (blk_y / p_o) + (left_over_y / 2);
	    int bottom_now = top_now + blk_x - 1;
	    int right_now = left_now + blk_y - 1;
	    /*printf("left_now: %i, right_now %i, top_now %i, bottom_now %i, p_o %i, blk_x %i, blk_y %i, width %i, height %i, Nlines_x %i, Nlines_y %i, left_over_x %i, left_over_y %i\n", left_now, right_now, top_now, bottom_now, p_o, blk_x, blk_y, w, h, Nlines_x, Nlines_y, left_over_x, left_over_y);*/
	    for (int i = top_now; i <= bottom_now; ++i) {
		for (int j = left_now; j <= right_now; ++j) {
		    number_t img1px = img1[i * img1_stride + j];
            	    number_t img2px = img2[i * img2_stride + j];
		    //printf("lala %f, %f\n", img1px, img2px);
		    motion_blocks_row[Nlines_y * (block_start_x_case - 1) + block_start_y_case - 1] += fabs(img1px - img2px);
		}
	    }
        }
    }

    for (int i = 0; i < height; ++i) {
        number_t accum_line = (number_t)0.0;
        for (int j = 0; j < width; ++j) {
            number_t img1px = img1[i * img1_stride + j];
            number_t img2px = img2[i * img2_stride + j];

            accum_line += fabs(img1px - img2px);

	    /*col_accum = MAX(floor(j / blk_y), 0);
            row_accum = MAX(floor(i / blk_x), 0);

            if ((col_accum < new_cols) && (row_accum < new_rows)) {
                motion_blocks_row[row_accum*new_cols + col_accum] +=  fabs(img1px - img2px);
            }*/

        }
                accum += accum_line;
    }

    /*printf("motion, blk_x %i, blk_y %i, width %i, height %i\n", blk_x, blk_y, width, height);*/

    printf("motion_blk: %d", frm_idx);
    for (int i = 0; i < Nlines_x * Nlines_y; i++) {
    	printf(" %f", motion_blocks_row[i] / (blk_x * blk_y));
    }
    printf("\n");
    fflush(stdout);

    return (number_t) (accum / (width * height));
}

/**
 * Note: ref_stride and dis_stride are in terms of bytes
 */
int compute_motion(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, int frm_idx, int blk_x, int blk_y, int new_rows, int new_cols, int p_o, int Nlines_x, int Nlines_y, int left_over_x, int left_over_y)
{

    if (ref_stride % sizeof(number_t) != 0)
    {
        printf("error: ref_stride %% sizeof(number_t) != 0, ref_stride = %d, sizeof(number_t) = %lu.\n", ref_stride, sizeof(number_t));
        fflush(stdout);
        goto fail;
    }
    if (dis_stride % sizeof(number_t) != 0)
    {
        printf("error: dis_stride %% sizeof(number_t) != 0, dis_stride = %d, sizeof(number_t) = %lu.\n", dis_stride, sizeof(number_t));
        fflush(stdout);
        goto fail;
    }
    // stride for vmaf_image_sad_c is in terms of (sizeof(number_t) bytes)
    *score = vmaf_image_sad_c(ref, dis, w, h, ref_stride / sizeof(number_t), dis_stride / sizeof(number_t), frm_idx, blk_x, blk_y, new_rows, new_cols, p_o, Nlines_x, Nlines_y, left_over_x, left_over_y);

    return 0;

fail:
    return 1;
}

int motion(const char *ref_path, int w, int h, const char *fmt, int blk_x_val, int blk_y_val, int p_o)
{
    double score = 0;
    FILE *ref_rfile = 0;
    number_t *ref_buf = 0;
    number_t *prev_blur_buf = 0;
    number_t *blur_buf = 0;
    number_t *temp_buf = 0;
    size_t data_sz;
    int stride;
    int ret = 1;
    /*int p_o = 2;*/
    int Nlines_x, Nlines_y, left_over_x, left_over_y;
    int blk_x, blk_y, new_cols, new_rows;
    int row_accum, col_accum;
    int i, ix_last, iy_last, ix_last_p_o, iy_last_p_o;

    if (w <= 0 || h <= 0 || (size_t)w > ALIGN_FLOOR(INT_MAX) / sizeof(number_t))
    {
        goto fail_or_end;
    }

    stride = ALIGN_CEIL(w * sizeof(number_t));

    if ((size_t)h > SIZE_MAX / stride)
    {
        goto fail_or_end;
    }

    data_sz = (size_t)stride * h;

    if (!(ref_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for ref_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }
    if (!(prev_blur_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for prev_blur_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }
    if (!(blur_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for blur_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }
    if (!(temp_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for temp_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    if (!(ref_rfile = fopen(ref_path, "rb")))
    {
        printf("error: fopen ref_path %s failed.\n", ref_path);
        fflush(stdout);
        goto fail_or_end;
    }

    size_t offset;
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv420p10le"))
    {
        if ((w * h) % 2 != 0)
        {
            printf("error: (w * h) %% 2 != 0, w = %d, h = %d.\n", w, h);
            fflush(stdout);
            goto fail_or_end;
        }
        offset = w * h / 2;
    }
    else if (!strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv422p10le"))
    {
        offset = w * h;
    }
    else if (!strcmp(fmt, "yuv444p") || !strcmp(fmt, "yuv444p10le"))
    {
        offset = w * h * 2;
    }
    else
    {
        printf("error: unknown format %s.\n", fmt);
        fflush(stdout);
        goto fail_or_end;
    }

    int frm_idx = 0;
    while (1)
    {
        // ref read y
        if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
        {
            ret = read_image_b(ref_rfile, ref_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
        }
        else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
        {
            ret = read_image_w(ref_rfile, ref_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
        }
        else
        {
            printf("error: unknown format %s.\n", fmt);
            fflush(stdout);
            goto fail_or_end;
        }
        if (ret)
        {
            if (feof(ref_rfile))
            {
                ret = 0; // OK if end of file
            }
            goto fail_or_end;
        }

        // filter
        // apply filtering (to eliminate effects film grain)
        // stride input to convolution_f32_c is in terms of (sizeof(number_t) bytes)
        // since stride = ALIGN_CEIL(w * sizeof(number_t)), stride divides sizeof(number_t)
        convolution_f32_c(FILTER_5, 5, ref_buf, blur_buf, temp_buf, w, h, stride / sizeof(number_t), stride / sizeof(number_t));

        // compute
        if (frm_idx == 0)
        {
            score = 0.0;
        }
        else
        {
            int left   = 0;
            int top    = 0;
            int right  = w;
            int bottom = h; 

            if ((blk_x_val == -1) || (blk_x_val == -1)) {
            	blk_x = bottom - top;
            	blk_y = right - left;
            }
            else {
            	blk_x = MIN(blk_x_val, bottom - top);
            	blk_y = MIN(blk_y_val, right - left);
            }

            i = top;
	    new_rows = 0;
            while (i <= bottom - blk_x) {
            	new_rows = new_rows + 1;
            	ix_last = i;
            	i = i + blk_x; 
            }

            i = left;
            new_cols = 0;
            while (i <= right - blk_y) {
            	new_cols = new_cols + 1;
            	iy_last = i;
            	i = i + blk_y; 
            }

            i = top;
            Nlines_x = 0;
            while (i <= bottom - blk_x) {
            	Nlines_x = Nlines_x + 1;
            	ix_last_p_o = i;
            	i = i + blk_x / p_o; 
            }

            i = left;
            Nlines_y = 0;
            while (i <= right - blk_y) {
            	Nlines_y = Nlines_y + 1;
            	iy_last_p_o = i;
            	i = i + blk_y / p_o; 
            }

            left_over_x = bottom - (MAX(ix_last, ix_last_p_o) + blk_x);
            left_over_y = right - (MAX(iy_last, iy_last_p_o) + blk_y);

            /*printf("new_rows: %i", new_rows);
            printf(" new_cols: %i", new_cols);
            printf(" Nlines_x: %i", Nlines_x);
            printf(" Nlines_y: %i\n", Nlines_y);*/

            /*new_cols = floor((right - left) / blk_y);
            new_rows = floor((bottom - top) / blk_x);
            Nlines_x = new_rows + (p_o - 1) * (new_rows - 1);
            Nlines_y = new_cols + (p_o - 1) * (new_cols - 1);*/
	
            if ((ret = compute_motion(prev_blur_buf, blur_buf, w, h, stride, stride, &score, frm_idx, blk_x, blk_y, new_rows, new_cols, p_o, Nlines_x, Nlines_y, left_over_x, left_over_y)))
            {
                printf("error: compute_motion failed.\n");
                fflush(stdout);
                goto fail_or_end;
            }
        }

        // copy to prev_buf
        memcpy(prev_blur_buf, blur_buf, data_sz);

        // print
        printf("motion: %d %f\n", frm_idx, score);
        fflush(stdout);

        /* remember to print for frame 0 and motion_blk as well */
        if (frm_idx == 0) {
		printf("motion_blk: %d", 0);
		for (int i = 0; i < Nlines_x * Nlines_y; i++) {
	       		printf(" %f", 0.0f);
		}
		printf("\n");
		fflush(stdout);
	}

        // ref skip u and v
        if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
        {
            if (fread(temp_buf, 1, offset, ref_rfile) != (size_t)offset)
            {
                printf("error: ref fread u and v failed.\n");
                fflush(stdout);
                goto fail_or_end;
            }
        }
        else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
        {
            if (fread(temp_buf, 2, offset, ref_rfile) != (size_t)offset)
            {
                printf("error: ref fread u and v failed.\n");
                fflush(stdout);
                goto fail_or_end;
            }
        }
        else
        {
            printf("error: unknown format %s.\n", fmt);
            fflush(stdout);
            goto fail_or_end;
        }

        frm_idx++;
    }

    ret = 0;

fail_or_end:
    if (ref_rfile)
    {
        fclose(ref_rfile);
    }
    aligned_free(ref_buf);
    aligned_free(prev_blur_buf);
    aligned_free(blur_buf);
    aligned_free(temp_buf);

    return ret;
}
