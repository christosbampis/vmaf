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

#include "common_christos/alloc.h"
#include "common_christos/file_io.h"
#include "common_christos/convolution.h"
#include "common_christos/convolution_internal.h"
#include "common_christos/tools.h"
#include "motion_tools.h"
#include "all_options.h"
#include "vif_options.h"
#include "adm_options.h"
#include "adm_tools.h"

#ifdef ALL_OPT_SINGLE_PRECISION
    typedef float number_t;

    #define read_image_b       read_image_b2s
    #define read_image_w       read_image_w2s
    #define convolution_f32_c  convolution_f32_c_s
    #define offset_image       offset_image_s
    #define FILTER_5           FILTER_5_s
    /*int compute_adm(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor);*/
    int compute_adm(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int Nblocks[Nscales], int px_stride[Nscales], int left[Nscales], int right[Nscales], int top[Nscales], int bottom[Nscales], int new_rows[Nscales], int new_cols[Nscales], int sum_of_Nblocks[Nscales], int frm_idx, int p_o, int Nlines_x[Nscales], int Nlines_y[Nscales], float adm_num_blk[Nscales], float adm_den_blk[Nscales], int left_over_x[Nscales], int left_over_y[Nscales], int is_it_diff);
    int compute_ansnr(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_psnr, double peak, double psnr_max);
    int compute_vif(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, int frm_idx, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int new_rows[Nscales], int new_cols[Nscales], int Nlines_x[Nscales], int Nlines_y[Nscales], int p_o, int left_over_x[Nscales], int left_over_y[Nscales], int is_it_diff);
    int compute_motion(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, int frm_idx, int blk_x_val, int blk_y_val, int new_rows, int new_cols, int p_o, int Nlines_x, int Nlines_y, int left_over_x, int left_over_y);

#else
    typedef double number_t;

    #define read_image_b       read_image_b2d
    #define read_image_w       read_image_w2d
    #define convolution_f32_c convolution_f32_c_d
    #define offset_image       offset_image_d
    #define FILTER_5           FILTER_5_d
    /*int compute_adm(const double *ref, const double *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor);*/
    int compute_adm(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int Nblocks[Nscales], int px_stride[Nscales], int left[Nscales], int right[Nscales], int top[Nscales], int bottom[Nscales], int new_rows[Nscales], int new_cols[Nscales], int sum_of_Nblocks[Nscales], int frm_idx, int p_o, int Nlines_x[Nscales], int Nlines_y[Nscales], float adm_num_blk[Nscales], float adm_den_blk[Nscales], int left_over_x[Nscales], int left_over_y[Nscales], int is_it_diff);
    int compute_ansnr(const double *ref, const double *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_psnr, double peak, double psnr_max);
    int compute_vif(const double *ref, const double *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, int frm_idx, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int new_rows[Nscales], int new_cols[Nscales], int Nlines_x[Nscales], int Nlines_y[Nscales], int p_o, int left_over_x[Nscales], int left_over_y[Nscales]);
    int compute_motion(const double *ref, const double *dis, int w, int h, int ref_stride, int dis_stride, double *score, int frm_idx, int blk_x_val, int blk_y_val, int new_rows, int new_cols, int p_o, int Nlines_x, int Nlines_y, int left_over_x, int left_over_y);

#endif

int all(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int start_x, int start_y, int roi_x, int roi_y, int blk_x_val, int blk_y_val, int p_o, int is_it_diff)
{
    double score = 0;
    double scores[4*2];
    double score_num = 0;
    double score_den = 0;
    double score_psnr = 0;
    number_t *ref_buf = 0;
    number_t *dis_buf = 0;

    number_t *ref_diff_buf = 0;
    number_t *dis_diff_buf = 0;

    number_t *prev_ref_buf = 0;
    number_t *prev_dis_buf = 0;
    number_t *temp_buf = 0;

    number_t *prev_blur_buf = 0;
    number_t *blur_buf = 0;

    FILE *ref_rfile = 0;
    FILE *dis_rfile = 0;
    size_t data_sz;
    int stride;
    int ret = 1;

    int Nscales = 4;
    /*int p_o = 2;*/
    int px_stride[Nscales], sum_of_Nblocks[Nscales];
    int left[Nscales], right[Nscales], top[Nscales], bottom[Nscales], Nblocks[Nscales];
    int left_adm[Nscales], right_adm[Nscales], top_adm[Nscales], bottom_adm[Nscales], Nblocks_adm[Nscales];
    int new_rows[Nscales], new_cols[Nscales], blk_x[Nscales], blk_y[Nscales];
    int new_rows_adm[Nscales], new_cols_adm[Nscales], blk_adm_x[Nscales], blk_adm_y[Nscales];
    int Nlines_x[Nscales], Nlines_y[Nscales], left_over_x[Nscales], left_over_y[Nscales], left_over_x_adm[Nscales], left_over_y_adm[Nscales];
    int Nlines_adm_x[Nscales], Nlines_adm_y[Nscales];
    float num_now, den_now;
    int i, ix_last, iy_last, ix_last_p_o, iy_last_p_o, ix_last_adm, iy_last_adm, ix_last_p_o_adm, iy_last_p_o_adm;

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
    if (!(dis_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for dis_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // prev_ref_buf
    if (!(prev_ref_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for prev_ref_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // ref_diff_buf
    if (!(ref_diff_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for ref_diff_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // dis_diff_buf
    if (!(dis_diff_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for dis_diff_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // prev_blur_buf
    if (!(prev_blur_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for prev_blur_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // blur_buf
    if (!(blur_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for blur_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // prev_dis_buf
    if (!(prev_dis_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for prev_dis_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    // use temp_buf for convolution_f32_c, and fread u and v
    if (!(temp_buf = aligned_malloc(data_sz * 2, MAX_ALIGN)))
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
    if (!(dis_rfile = fopen(dis_path, "rb")))
    {
        printf("error: fopen dis_path %s failed.\n", dis_path);
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
        // read ref y
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

        // read dis y
        if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
        {
            ret = read_image_b(dis_rfile, dis_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
        }
        else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
        {
            ret = read_image_w(dis_rfile, dis_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
        }
        else
        {
            printf("error: unknown format %s.\n", fmt);
            fflush(stdout);
            goto fail_or_end;
        }
        if (ret)
        {
            if (feof(dis_rfile))
            {
                ret = 0; // OK if end of file
            }
            goto fail_or_end;
        }

	int stride_block = ALIGN_CEIL(roi_y * sizeof(number_t));
	int w_now = w;
	int h_now = h;

	int w_adm_now = w;
	int h_adm_now = h;

	int buf_stride = ALIGN_CEIL(((w + 1) / 2) * sizeof(number_t));

	for (int scale_ind = 0; scale_ind < Nscales; scale_ind++) { sum_of_Nblocks[scale_ind] = 0; };

	/* remember that adm starts of at the lower scale */
	int blk_x_val_adm = (blk_x_val + 1) / 2;
	int blk_y_val_adm = (blk_y_val + 1) / 2;

	for (int scale_ind = 0; scale_ind < Nscales; scale_ind++) { 

		w_adm_now = (w_adm_now + 1) / 2;
		h_adm_now = (h_adm_now + 1) / 2;

		Nblocks[scale_ind] = 0; 
		px_stride[scale_ind] = buf_stride / sizeof(number_t);
		/*left[scale_ind]   = w_now * ADM_BORDER_FACTOR - 0.5;
	    	top[scale_ind]    = h_now * ADM_BORDER_FACTOR - 0.5;*/
		left[scale_ind]   = MAX(w_now * ADM_BORDER_FACTOR - 0.5, 0); 
	    	top[scale_ind]    = MAX(h_now * ADM_BORDER_FACTOR - 0.5, 0);

		left_adm[scale_ind]   = MAX(w_adm_now * ADM_BORDER_FACTOR - 0.5, 0); 
	    	top_adm[scale_ind]    = MAX(h_adm_now * ADM_BORDER_FACTOR - 0.5, 0);

		right[scale_ind]  = w_now - left[scale_ind];
		bottom[scale_ind] = h_now - top[scale_ind]; 

		right_adm[scale_ind]  = w_adm_now - left_adm[scale_ind];
		bottom_adm[scale_ind] = h_adm_now - top_adm[scale_ind]; 

		if ((blk_x_val == -1) || (blk_x_val == -1)) {
			blk_x[scale_ind] = bottom[scale_ind] - top[scale_ind];
			blk_y[scale_ind] = right[scale_ind] - left[scale_ind];
			blk_adm_x[scale_ind] = bottom_adm[scale_ind] - top_adm[scale_ind];
			blk_adm_y[scale_ind] = right_adm[scale_ind] - left_adm[scale_ind];
    		}
		else {
			if (scale_ind == 0) {
				blk_x[scale_ind] = MIN(blk_x_val, bottom[scale_ind] - top[scale_ind]);
				blk_y[scale_ind] = MIN(blk_y_val, right[scale_ind] - left[scale_ind]);

				blk_adm_x[scale_ind] = MIN(blk_x_val_adm, bottom_adm[scale_ind] - top_adm[scale_ind]);
				blk_adm_y[scale_ind] = MIN(blk_y_val_adm, right_adm[scale_ind] - left_adm[scale_ind]);
			}
			else {
				blk_x[scale_ind] = MIN((blk_x[scale_ind-1] + 1) / 2, bottom[scale_ind] - top[scale_ind]);
				blk_y[scale_ind] = MIN((blk_y[scale_ind-1] + 1) / 2, right[scale_ind] - left[scale_ind]);

				blk_adm_x[scale_ind] = MIN((blk_adm_x[scale_ind-1] + 1) / 2, bottom_adm[scale_ind] - top_adm[scale_ind]);
				blk_adm_y[scale_ind] = MIN((blk_adm_y[scale_ind-1] + 1) / 2, right_adm[scale_ind] - left_adm[scale_ind]);
			}
		}

		i = top[scale_ind];
		new_rows[scale_ind] = 0;
		while (i <= bottom[scale_ind] - blk_x[scale_ind]) {
			new_rows[scale_ind] = new_rows[scale_ind] + 1;
			ix_last = i;
			i = i + blk_x[scale_ind]; 
		}

		i = left[scale_ind];
		new_cols[scale_ind] = 0;
		while (i <= right[scale_ind] - blk_y[scale_ind]) {
			new_cols[scale_ind] = new_cols[scale_ind] + 1;
			iy_last = i;
			i = i + blk_y[scale_ind]; 
		}

		i = top[scale_ind];
		Nlines_x[scale_ind] = 0;
		while (i <= bottom[scale_ind] - blk_x[scale_ind]) {
			Nlines_x[scale_ind] = Nlines_x[scale_ind] + 1;
			ix_last_p_o = i;
			i = i + blk_x[scale_ind] / p_o; 
		}

		i = left[scale_ind];
		Nlines_y[scale_ind] = 0;
		while (i <= right[scale_ind] - blk_y[scale_ind]) {
			Nlines_y[scale_ind] = Nlines_y[scale_ind] + 1;
			iy_last_p_o = i;
			i = i + blk_y[scale_ind] / p_o; 
		}

		left_over_x[scale_ind] = bottom[scale_ind] - (MAX(ix_last, ix_last_p_o) + blk_x[scale_ind]);
		left_over_y[scale_ind] = right[scale_ind] - (MAX(iy_last, iy_last_p_o) + blk_y[scale_ind]);

		i = top_adm[scale_ind];
		new_rows_adm[scale_ind] = 0;
		while (i <= bottom_adm[scale_ind] - blk_adm_x[scale_ind]) {
			new_rows_adm[scale_ind] = new_rows_adm[scale_ind] + 1;
			ix_last_adm = i;
			i = i + blk_adm_x[scale_ind]; 
		}

		i = left_adm[scale_ind];
		new_cols_adm[scale_ind] = 0;
		while (i <= right_adm[scale_ind] - blk_adm_y[scale_ind]) {
			new_cols_adm[scale_ind] = new_cols_adm[scale_ind] + 1;
			iy_last_adm = i;
			i = i + blk_adm_y[scale_ind]; 
		}

		i = top_adm[scale_ind];
		Nlines_adm_x[scale_ind] = 0;
		while (i <= bottom_adm[scale_ind] - blk_adm_x[scale_ind]) {
			Nlines_adm_x[scale_ind] = Nlines_adm_x[scale_ind] + 1;
			ix_last_p_o_adm = i;
			i = i + blk_adm_x[scale_ind] / p_o; 
		}

		i = left_adm[scale_ind];
		Nlines_adm_y[scale_ind] = 0;
		while (i <= right_adm[scale_ind] - blk_adm_y[scale_ind]) {
			Nlines_adm_y[scale_ind] = Nlines_adm_y[scale_ind] + 1;
			iy_last_p_o_adm = i;
			i = i + blk_adm_y[scale_ind] / p_o; 
		}

		left_over_x_adm[scale_ind] = bottom_adm[scale_ind] - (MAX(ix_last_adm, ix_last_p_o_adm) + blk_adm_x[scale_ind]);
		left_over_y_adm[scale_ind] = right_adm[scale_ind] - (MAX(iy_last_adm, iy_last_p_o_adm) + blk_adm_y[scale_ind]);

	    	/*new_cols[scale_ind] = floor((right[scale_ind] - left[scale_ind]) / blk_y[scale_ind]);
	    	new_rows[scale_ind] = floor((bottom[scale_ind] - top[scale_ind]) / blk_x[scale_ind]);

	    	new_cols_adm[scale_ind] = floor((right_adm[scale_ind] - left_adm[scale_ind]) / blk_adm_y[scale_ind]);
	    	new_rows_adm[scale_ind] = floor((bottom_adm[scale_ind] - top_adm[scale_ind]) / blk_adm_x[scale_ind]);

		Nlines_x[scale_ind] = new_rows[scale_ind] + (p_o - 1) * (new_rows[scale_ind] - 1);
		Nlines_y[scale_ind] = new_cols[scale_ind] + (p_o - 1) * (new_cols[scale_ind] - 1);

		Nlines_adm_x[scale_ind] = new_rows_adm[scale_ind] + (p_o - 1) * (new_rows_adm[scale_ind] - 1);
		Nlines_adm_y[scale_ind] = new_cols_adm[scale_ind] + (p_o - 1) * (new_cols_adm[scale_ind] - 1);*/

	    	/*printf("blk_x %i, blk_y %i, new_cols %i, new_rows %i\n", blk_x[scale_ind], blk_y[scale_ind], new_cols[scale_ind], new_rows[scale_ind]);*/
		Nblocks[scale_ind] = Nlines_x[scale_ind] * Nlines_y[scale_ind];
		Nblocks_adm[scale_ind] = Nlines_adm_x[scale_ind] * Nlines_adm_y[scale_ind];

		w_now = (w_now + 1) / 2;
		h_now = (h_now + 1) / 2;

	};

	float adm_num_blk[Nblocks_adm[0]], adm_den_blk[Nblocks_adm[0]];
    	for (int i = 0; i < Nblocks[0]; ++i) {
		adm_num_blk[i] = 0.0;
		adm_den_blk[i] = 0.0;
    	}


	if (is_it_diff == 1) {
		if (frm_idx > 0) {
			frame_differ(ref_buf, prev_ref_buf, ref_diff_buf, w, h, stride / sizeof(number_t));
			frame_differ(dis_buf, prev_dis_buf, dis_diff_buf, w, h, stride / sizeof(number_t));	
		}
	}

	memcpy(prev_ref_buf, ref_buf, data_sz);
	memcpy(prev_dis_buf, dis_buf, data_sz);

	if (is_it_diff == 1) {

		if (frm_idx == 0)
		{
		    score = 0.0;
		    score_num = 0.0;
		    score_den = 0.0;
		    for(int scale=0; scale<Nscales; scale++){
		    	scores[2*scale] = 0.0;
		    	scores[2*scale+1] = 0.0;
		    }
		}
		else
		{
			/*printf("/////////\n");
			for (int i = 0; i < 10; ++i) {
				for (int j = 0; j < 10; ++j) {
					number_t img_px_after = ref_buf[i * stride / sizeof(number_t) + j];
					number_t img_px_before = prev_ref_buf[i * stride / sizeof(number_t) + j];
					number_t img_px_diff = ref_diff_buf[i * stride / sizeof(number_t) + j];
					printf("i: %f, i+1: %f, diff: %f\n", img_px_after, img_px_before, img_px_diff);
				}
			}*/

			/* =========== adm ============== */
			if ((ret = compute_adm(ref_diff_buf, dis_diff_buf, w, h, stride_block, stride_block, &score, &score_num, &score_den, scores, ADM_BORDER_FACTOR, Nscales, blk_adm_x, blk_adm_y, Nblocks_adm, px_stride, left_adm, right_adm, top_adm, bottom_adm, new_rows_adm, new_cols_adm, sum_of_Nblocks, frm_idx, p_o, Nlines_x, Nlines_y, adm_num_blk, adm_den_blk, left_over_x_adm, left_over_y_adm, is_it_diff)))
			{
				printf("error: compute_adm_diff failed.\n");
				fflush(stdout);
				goto fail_or_end;
			}

		}
	}

	else {
		/* =========== adm ============== */
		if ((ret = compute_adm(ref_buf, dis_buf, w, h, stride_block, stride_block, &score, &score_num, &score_den, scores, ADM_BORDER_FACTOR, Nscales, blk_adm_x, blk_adm_y, Nblocks_adm, px_stride, left_adm, right_adm, top_adm, bottom_adm, new_rows_adm, new_cols_adm, sum_of_Nblocks, frm_idx, p_o, Nlines_x, Nlines_y, adm_num_blk, adm_den_blk, left_over_x_adm, left_over_y_adm, is_it_diff)))
		{
			printf("error: compute_adm failed.\n");
			fflush(stdout);
			goto fail_or_end;
		}
	}

	if (is_it_diff == 1) {
		printf("adm_diff: %d %f\n", frm_idx, score);
	}
	else {
		printf("adm: %d %f\n", frm_idx, score);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
		printf("adm_num_diff: %d %f\n", frm_idx, score_num);
	}
	else {
		printf("adm_num: %d %f\n", frm_idx, score_num);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
		printf("adm_den_diff: %d %f\n", frm_idx, score_den);
	}
	else {
		printf("adm_den: %d %f\n", frm_idx, score_den);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
		for(int scale=0; scale<Nscales; scale++){
			printf("adm_num_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale]);
			printf("adm_den_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}
	else {
		for(int scale=0; scale<Nscales; scale++){
			printf("adm_num_scale%d: %d %f\n", scale, frm_idx, scores[2*scale]);
			printf("adm_den_scale%d: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}
	fflush(stdout);

        /* =========== ansnr ============== */
        /*if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
        {
            // max psnr 60.0 for 8-bit per Ioannis
            ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.0, 60.0);
        }
        else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
        {
            // 10 bit gets normalized to 8 bit, peak is 1023 / 4.0 = 255.75
            // max psnr 72.0 for 10-bit per Ioannis
            ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.75, 72.0);
        }
        else
        {
            printf("error: unknown format %s.\n", fmt);
            fflush(stdout);
            goto fail_or_end;
        }
        if (ret)
        {
            printf("error: compute_ansnr failed.\n");
            fflush(stdout);
            goto fail_or_end;
        }

        printf("ansnr: %d %f\n", frm_idx, score);
        fflush(stdout);
        printf("anpsnr: %d %f\n", frm_idx, score_psnr);
        fflush(stdout);*/

        /* =========== motion ============== */

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

           if ((ret = compute_motion(prev_blur_buf, blur_buf, w, h, stride, stride, &score, frm_idx, blk_x[0], blk_y[0], new_rows[0], new_cols[0], p_o, Nlines_x[0], Nlines_y[0], left_over_x[0], left_over_y[0])))
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

        /* =========== vif ============== */
        // compute vif last, because its input ref/dis might be offset by -128
	/*for (int i = 0; i < 4; i++) {
		printf("new cols:%i, new rows: %i\n", new_cols[i], new_rows[i]);
	}*/

	/* have to change this if I want to use is_it_diff for both instances */
	
	if (is_it_diff == 1) {

		if (frm_idx == 0)
		{
		    score = 0.0;
		    score_num = 0.0;
		    score_den = 0.0;
		    for(int scale=0; scale<Nscales; scale++){
		    	scores[2*scale] = 0.0;
		    	scores[2*scale+1] = 0.0;
		    }
		}
		else
		{

			if ((ret = compute_vif(ref_diff_buf, dis_diff_buf, w, h, stride, stride, &score, &score_num, &score_den, scores, frm_idx, Nscales, blk_x, blk_y, new_rows, new_cols, Nlines_x, Nlines_y, p_o, left_over_x, left_over_y, is_it_diff)))
			{
				printf("error: compute_vif_diff failed.\n");
				fflush(stdout);
				goto fail_or_end;
			}
		}

	}

	else {

		if ((ret = compute_vif(ref_buf, dis_buf, w, h, stride, stride, &score, &score_num, &score_den, scores, frm_idx, Nscales, blk_x, blk_y, new_rows, new_cols, Nlines_x, Nlines_y, p_o, left_over_x, left_over_y, is_it_diff)))
		{
			printf("error: compute_vif failed.\n");
			fflush(stdout);
			goto fail_or_end;
		}

	}

	if (is_it_diff == 1) {
		printf("vif_diff: %d %f\n", frm_idx, score);
	}
	else {
		printf("vif: %d %f\n", frm_idx, score);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
	printf("vif_num_diff: %d %f\n", frm_idx, score_num);
	}
	else {
		printf("vif_num: %d %f\n", frm_idx, score_num);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
	printf("vif_den_diff: %d %f\n", frm_idx, score_den);
	}
	else {
		printf("vif_den: %d %f\n", frm_idx, score_den);
	}
	fflush(stdout);

	if (is_it_diff == 1) {
		for(int scale=0;scale<4;scale++){
			printf("vif_num_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale]);
			printf("vif_den_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}
	else {
		for(int scale=0;scale<4;scale++){
			printf("vif_num_scale%d: %d %f\n", scale, frm_idx, scores[2*scale]);
			printf("vif_den_scale%d: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
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

        // dis skip u and v
        if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
        {
            if (fread(temp_buf, 1, offset, dis_rfile) != (size_t)offset)
            {
                printf("error: dis fread u and v failed.\n");
                fflush(stdout);
                goto fail_or_end;
            }
        }
        else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
        {
            if (fread(temp_buf, 2, offset, dis_rfile) != (size_t)offset)
            {
                printf("error: dis fread u and v failed.\n");
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

	/* remember to print for frame 0 as well */
        if (is_it_diff == 1 && frm_idx == 0) {

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("vif_num_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("vif_den_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("vif_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

		printf("motion_blk: %i", frm_idx);
		for (int i = 0; i < Nlines_x[0] * Nlines_y[0]; i++) {
		       	printf(" %f", 0.0f);
		}
		printf("\n");
		fflush(stdout);

		printf("adm_num_diff_blk: %i", frm_idx);
		for (int i = 0; i < Nlines_x[0] * Nlines_y[0]; i++) {
		       	printf(" %f", 0.0f);
		}
		printf("\n");
		fflush(stdout);

		printf("adm_den_diff_blk: %i", frm_idx);
		for (int i = 0; i < Nlines_x[0] * Nlines_y[0]; i++) {
		       	printf(" %f", 0.0f);
		}
		printf("\n");
		fflush(stdout);

		printf("adm_diff_blk: %i", frm_idx);
		for (int i = 0; i < Nlines_x[0] * Nlines_y[0]; i++) {
		       	printf(" %f", 0.0f);
		}
		printf("\n");
		fflush(stdout);

		printf("adm_sm_diff_blk: %i %f\n", frm_idx, score);
		fflush(stdout);

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("adm_num_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("adm_den_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

		for(int scale_ind = 0; scale_ind < Nscales; scale_ind++){
			printf("adm_scale%i_diff_blk: %i", scale_ind, frm_idx);
			for (int i = 0; i < Nlines_x[scale_ind] * Nlines_y[scale_ind]; i++) {
		       		printf(" %f", 0.0f);
			}
			printf("\n");
		}
		fflush(stdout);

	}

        frm_idx++;
    }

    ret = 0;

fail_or_end:
    if (ref_rfile)
    {
        fclose(ref_rfile);
    }
    if (dis_rfile)
    {
        fclose(dis_rfile);
    }
    aligned_free(ref_buf);
    aligned_free(dis_buf);

    aligned_free(ref_diff_buf);
    aligned_free(dis_diff_buf);

    aligned_free(prev_blur_buf);
    aligned_free(blur_buf);

    aligned_free(prev_ref_buf);
    aligned_free(prev_dis_buf);
    aligned_free(temp_buf);

    return ret;
}
