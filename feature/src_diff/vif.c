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
#include "vif_options.h"
#include "vif_tools.h"
#include "common_christos/tools.h"

#ifdef VIF_OPT_SINGLE_PRECISION
  typedef float number_t;

  #define read_image_b       read_image_b2s
  #define read_image_w       read_image_w2s
  #define vif_filter1d_table vif_filter1d_table_s
  #define vif_filter1d       vif_filter1d_s
  #define vif_filter2d_table vif_filter2d_table_s
  #define vif_filter2d       vif_filter2d_s
  #define vif_dec2           vif_dec2_s
  #define vif_sum_blk        vif_sum_s_blk
  #define vif_sum            vif_sum_s
  #define vif_xx_yy_xy       vif_xx_yy_xy_s
  #define vif_statistic      vif_statistic_s
#else
  typedef double number_t;

  #define read_image_b       read_image_b2d
  #define read_image_w       read_image_w2d
  #define vif_filter1d_table vif_filter1d_table_d
  #define vif_filter1d       vif_filter1d_d
  #define vif_filter2d_table vif_filter2d_table_d
  #define vif_filter2d       vif_filter2d_d
  #define vif_dec2           vif_dec2_d
  #define vif_sum            vif_sum_d
  #define vif_xx_yy_xy       vif_xx_yy_xy_d
  #define vif_statistic      vif_statistic_d
#endif

int compute_vif(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, int frm_idx, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int new_rows[Nscales], int new_cols[Nscales], int Nlines_x[Nscales], int Nlines_y[Nscales], int p_o, int left_over_x[Nscales], int left_over_y[Nscales], int is_it_diff)
{

    number_t *data_buf = 0;
    char *data_top;

    number_t *ref_scale;
    number_t *dis_scale;
    number_t *ref_sq;
    number_t *dis_sq;
    number_t *ref_dis;

    number_t *mu1;
    number_t *mu2;
    number_t *mu1_sq;
    number_t *mu2_sq;
    number_t *mu1_mu2;
    number_t *ref_sq_filt;
    number_t *dis_sq_filt;
    number_t *ref_dis_filt;
    number_t *num_array;
    number_t *den_array;
    number_t *tmpbuf;

    /* Offset pointers to adjust for convolution border handling. */
    number_t *mu1_adj = 0;
    number_t *mu2_adj = 0;

#ifdef VIF_OPT_DEBUG_DUMP
    number_t *mu1_sq_adj;
    number_t *mu2_sq_adj;
    number_t *mu1_mu2_adj;
    number_t *ref_sq_filt_adj;
    number_t *dis_sq_filt_adj;
    number_t *ref_dis_filt_adj = 0;
#endif

    number_t *num_array_adj = 0;
    number_t *den_array_adj = 0;

    /* Special handling of first scale. */
    const number_t *curr_ref_scale = ref;
    const number_t *curr_dis_scale = dis;
    int curr_ref_stride = ref_stride;
    int curr_dis_stride = dis_stride;

    int buf_stride = ALIGN_CEIL(w * sizeof(number_t));
    size_t buf_sz_one = (size_t)buf_stride * h;

    double num = 0;
    double den = 0;

    int scale;
    int ret = 1;

    if (SIZE_MAX / buf_sz_one < 15)
    {
        printf("error: SIZE_MAX / buf_sz_one < 15, buf_sz_one = %lu.\n", buf_sz_one);
        fflush(stdout);
        goto fail_or_end;
    }

    if (!(data_buf = aligned_malloc(buf_sz_one * 16, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for data_buf.\n");
        fflush(stdout);
        goto fail_or_end;
    }

    data_top = (char *)data_buf;

    ref_scale = (number_t *)data_top; data_top += buf_sz_one;
    dis_scale = (number_t *)data_top; data_top += buf_sz_one;
    ref_sq    = (number_t *)data_top; data_top += buf_sz_one;
    dis_sq    = (number_t *)data_top; data_top += buf_sz_one;
    ref_dis   = (number_t *)data_top; data_top += buf_sz_one;
    mu1          = (number_t *)data_top; data_top += buf_sz_one;
    mu2          = (number_t *)data_top; data_top += buf_sz_one;
    mu1_sq       = (number_t *)data_top; data_top += buf_sz_one;
    mu2_sq       = (number_t *)data_top; data_top += buf_sz_one;
    mu1_mu2      = (number_t *)data_top; data_top += buf_sz_one;
    ref_sq_filt  = (number_t *)data_top; data_top += buf_sz_one;
    dis_sq_filt  = (number_t *)data_top; data_top += buf_sz_one;
    ref_dis_filt = (number_t *)data_top; data_top += buf_sz_one;
    num_array    = (number_t *)data_top; data_top += buf_sz_one;
    den_array    = (number_t *)data_top; data_top += buf_sz_one;
    tmpbuf    = (number_t *)data_top; data_top += buf_sz_one;

    float num_blk_sm = 0.0f;
    float den_blk_sm = 0.0f;

    for (scale = 0; scale < Nscales; ++scale)
    {
#ifdef VIF_OPT_DEBUG_DUMP
        char pathbuf[256];
#endif

#ifdef VIF_OPT_FILTER_1D
        const number_t *filter = vif_filter1d_table[scale];
        int filter_width       = vif_filter1d_width[scale];
#else
        const number_t *filter = vif_filter2d_table[scale];
        int filter_width       = vif_filter2d_width[scale];
#endif

#ifdef VIF_OPT_HANDLE_BORDERS
        int buf_valid_w = w;
        int buf_valid_h = h;

  #define ADJUST(x) x
#else
        int filter_adj  = filter_width / 2;
        int buf_valid_w = w - filter_adj * 2;
        int buf_valid_h = h - filter_adj * 2;

  #define ADJUST(x) ((number_t *)((char *)(x) + filter_adj * buf_stride + filter_adj * sizeof(number_t)))
#endif

        if (scale > 0)
        {
#ifdef VIF_OPT_FILTER_1D
            vif_filter1d(filter, curr_ref_scale, mu1, tmpbuf, w, h, curr_ref_stride, buf_stride, filter_width);
            vif_filter1d(filter, curr_dis_scale, mu2, tmpbuf, w, h, curr_dis_stride, buf_stride, filter_width);
#else
            vif_filter2d(filter, curr_ref_scale, mu1, w, h, curr_ref_stride, buf_stride, filter_width);
            vif_filter2d(filter, curr_dis_scale, mu2, w, h, curr_dis_stride, buf_stride, filter_width);
#endif
            mu1_adj = ADJUST(mu1);
            mu2_adj = ADJUST(mu2);

            vif_dec2(mu1_adj, ref_scale, buf_valid_w, buf_valid_h, buf_stride, buf_stride);
            vif_dec2(mu2_adj, dis_scale, buf_valid_w, buf_valid_h, buf_stride, buf_stride);

            w  = buf_valid_w / 2;
            h  = buf_valid_h / 2;
#ifdef VIF_OPT_HANDLE_BORDERS
            buf_valid_w = w;
            buf_valid_h = h;
#else
            buf_valid_w = w - filter_adj * 2;
            buf_valid_h = h - filter_adj * 2;
#endif
            curr_ref_scale = ref_scale;
            curr_dis_scale = dis_scale;

            curr_ref_stride = buf_stride;
            curr_dis_stride = buf_stride;
        }

#ifdef VIF_OPT_FILTER_1D
        vif_filter1d(filter, curr_ref_scale, mu1, tmpbuf, w, h, curr_ref_stride, buf_stride, filter_width);
        vif_filter1d(filter, curr_dis_scale, mu2, tmpbuf, w, h, curr_dis_stride, buf_stride, filter_width);
#else
        vif_filter2d(filter, curr_ref_scale, mu1, w, h, curr_ref_stride, buf_stride, filter_width);
        vif_filter2d(filter, curr_dis_scale, mu2, w, h, curr_dis_stride, buf_stride, filter_width);
#endif
        vif_xx_yy_xy(mu1, mu2, mu1_sq, mu2_sq, mu1_mu2, w, h, buf_stride, buf_stride, buf_stride, buf_stride, buf_stride);

        vif_xx_yy_xy(curr_ref_scale, curr_dis_scale, ref_sq, dis_sq, ref_dis, w, h, curr_ref_stride, curr_dis_stride, buf_stride, buf_stride, buf_stride);
#ifdef VIF_OPT_FILTER_1D
        vif_filter1d(filter, ref_sq, ref_sq_filt, tmpbuf, w, h, buf_stride, buf_stride, filter_width);
        vif_filter1d(filter, dis_sq, dis_sq_filt, tmpbuf, w, h, buf_stride, buf_stride, filter_width);
        vif_filter1d(filter, ref_dis, ref_dis_filt, tmpbuf, w, h, buf_stride, buf_stride, filter_width);
#else
        vif_filter2d(filter, ref_sq, ref_sq_filt, w, h, buf_stride, buf_stride, filter_width);
        vif_filter2d(filter, dis_sq, dis_sq_filt, w, h, buf_stride, buf_stride, filter_width);
        vif_filter2d(filter, ref_dis, ref_dis_filt, w, h, buf_stride, buf_stride, filter_width);
#endif
        vif_statistic(mu1_sq, mu2_sq, mu1_mu2, ref_sq_filt, dis_sq_filt, ref_dis_filt, num_array, den_array,
                      w, h, buf_stride, buf_stride, buf_stride, buf_stride, buf_stride, buf_stride, buf_stride, buf_stride);

        mu1_adj = ADJUST(mu1);
        mu2_adj = ADJUST(mu2);

#ifdef VIF_OPT_DEBUG_DUMP
        mu1_sq_adj  = ADJUST(mu1_sq);
        mu2_sq_adj  = ADJUST(mu2_sq);
        mu1_mu2_adj = ADJUST(mu1_mu2);

        ref_sq_filt_adj  = ADJUST(ref_sq_filt);
        dis_sq_filt_adj  = ADJUST(dis_sq_filt);
        ref_dis_filt_adj = ADJUST(ref_dis_filt);
#endif

        num_array_adj = ADJUST(num_array);
        den_array_adj = ADJUST(den_array);
#undef ADJUST

#ifdef VIF_OPT_DEBUG_DUMP
        sprintf(pathbuf, "stage/ref[%d].bin", scale);
        write_image(pathbuf, curr_ref_scale, w, h, curr_ref_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis[%d].bin", scale);
        write_image(pathbuf, curr_dis_scale, w, h, curr_dis_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mu1[%d].bin", scale);
        write_image(pathbuf, mu1_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mu2[%d].bin", scale);
        write_image(pathbuf, mu2_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mu1_sq[%d].bin", scale);
        write_image(pathbuf, mu1_sq_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mu2_sq[%d].bin", scale);
        write_image(pathbuf, mu2_sq_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mu1_mu2[%d].bin", scale);
        write_image(pathbuf, mu1_mu2_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/ref_sq_filt[%d].bin", scale);
        write_image(pathbuf, ref_sq_filt_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis_sq_filt[%d].bin", scale);
        write_image(pathbuf, dis_sq_filt_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/ref_dis_filt[%d].bin", scale);
        write_image(pathbuf, ref_dis_filt_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/num_array[%d].bin", scale);
        write_image(pathbuf, num_array_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/den_array[%d].bin", scale);
        write_image(pathbuf, den_array_adj, buf_valid_w, buf_valid_h, buf_stride, sizeof(number_t));
#endif

	float num_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];
	float den_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];
	float vif_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) { 
		num_scale_blocks_row[i] = 0.0f;
		den_scale_blocks_row[i] = 0.0f;
		vif_scale_blocks_row[i] = 0.0f;
	}

	vif_sum_blk(num_array_adj, buf_valid_w, buf_valid_h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], num_scale_blocks_row, Nlines_x[scale], Nlines_y[scale], p_o, left_over_x[scale], left_over_y[scale]);
	vif_sum_blk(den_array_adj, buf_valid_w, buf_valid_h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], den_scale_blocks_row, Nlines_x[scale], Nlines_y[scale], p_o, left_over_x[scale], left_over_y[scale]);

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) {
		if (den_scale_blocks_row[i] == 0.0) {
			vif_scale_blocks_row[i] = 1.0f;
		}
		else {		
			vif_scale_blocks_row[i] = num_scale_blocks_row[i] / den_scale_blocks_row[i];
		}
	}

	if (is_it_diff == 1) {
		printf("vif_num_scale%i_diff_blk: %i", scale, frm_idx);
	}
	else {
		printf("vif_num_scale%i_blk: %i", scale, frm_idx);
	}

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) {
		printf(" %f", num_scale_blocks_row[i]);
	}
	printf("\n");
	fflush(stdout);

	if (is_it_diff == 1) {
		printf("vif_den_scale%i_diff_blk: %i", scale, frm_idx);
	}
	else {
		printf("vif_den_scale%i_blk: %i", scale, frm_idx);
	}

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) {
		printf(" %f", den_scale_blocks_row[i]);
	}
	printf("\n");
	fflush(stdout);

	if (is_it_diff == 1) {
		printf("vif_scale%i_diff_blk: %i", scale, frm_idx);
	}
	else {
		printf("vif_scale%i_blk: %i", scale, frm_idx);
	}

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) {
		printf(" %f", vif_scale_blocks_row[i]);
	}
	printf("\n");
	fflush(stdout);

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) {

		num_blk_sm += num_scale_blocks_row[i];
		den_blk_sm += den_scale_blocks_row[i];

	}

        num = vif_sum(num_array_adj, buf_valid_w, buf_valid_h, buf_stride);
        den = vif_sum(den_array_adj, buf_valid_w, buf_valid_h, buf_stride);

        scores[2*scale] = num;
        scores[2*scale+1] = den;

#ifdef VIF_OPT_DEBUG_DUMP
        printf("num[%d]: %e\n", scale, num);
        printf("den[%d]: %e\n", scale, den);
#endif
    }

    if (is_it_diff == 1) {
    	printf("vif_sm_diff_blk: %i %f\n", frm_idx, num_blk_sm / den_blk_sm);
    }
    else {
 	printf("vif_sm_blk: %i %f\n", frm_idx, num_blk_sm / den_blk_sm);
    }
    fflush(stdout);

    *score_num = 0.0;
    *score_den = 0.0;
    for (scale = 0; scale < 4; ++scale)
    {
        *score_num += scores[2*scale];
        *score_den += scores[2*scale+1];
    }
    if (*score_den == 0.0)
    {
        *score = 1.0f;
    }
    else
    {
        *score = (*score_num) / (*score_den);
    }

    ret = 0;
fail_or_end:
    aligned_free(data_buf);
    return ret;
}

int vif(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int blk_x_val, int blk_y_val, int p_o, int is_it_diff)
{
    double score = 0;
    double scores[4*2];
    double score_num = 0;
    double score_den = 0;

    number_t *ref_buf = 0;
    number_t *dis_buf = 0;
    number_t *temp_buf = 0;

    number_t *ref_diff_buf = 0;
    number_t *dis_diff_buf = 0;

    number_t *prev_ref_buf = 0;
    number_t *prev_dis_buf = 0;

    FILE *ref_rfile = 0;
    FILE *dis_rfile = 0;
    size_t data_sz;
    int stride;
    int ret = 1;

    int i, ix_last, iy_last, ix_last_p_o, iy_last_p_o;
    int Nscales = 4;
    /*int p_o = 2;*/
    int left[Nscales], right[Nscales], top[Nscales], bottom[Nscales], Nblocks[Nscales];	
    int new_rows[Nscales], new_cols[Nscales], blk_x[Nscales], blk_y[Nscales];
    int Nlines_x[Nscales], Nlines_y[Nscales], left_over_x[Nscales], left_over_y[Nscales];

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

    // prev_dis_buf
    if (!(prev_dis_buf = aligned_malloc(data_sz, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for prev_dis_buf.\n");
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

	int w_now = w;
	int h_now = h;

	for (int scale_ind = 0; scale_ind < Nscales; scale_ind++) { 

		Nblocks[scale_ind] = 0; 
		/* could have defined this in terms of ADM (using ADM_BORDER_FACTOR) as done in all.c */
		left[scale_ind]   = 0; 
	    	top[scale_ind]    = 0;
		right[scale_ind]  = w_now - left[scale_ind];
		bottom[scale_ind] = h_now - top[scale_ind];  
		if ((blk_x_val == -1) || (blk_x_val == -1)) {
			blk_x[scale_ind] = bottom[scale_ind] - top[scale_ind];
			blk_y[scale_ind] = right[scale_ind] - left[scale_ind];
    		}
		else {
			if (scale_ind == 0) {
				blk_x[scale_ind] = MIN(blk_x_val, bottom[scale_ind] - top[scale_ind]);
				blk_y[scale_ind] = MIN(blk_y_val, right[scale_ind] - left[scale_ind]);
			}
			else {
				blk_x[scale_ind] = MIN((blk_x[scale_ind-1] + 1) / 2, bottom[scale_ind] - top[scale_ind]);
				blk_y[scale_ind] = MIN((blk_y[scale_ind-1] + 1) / 2, right[scale_ind] - left[scale_ind]);
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
		/*printf("iy_last: %i, iy_last_p_o %i mx %i loy %i\n", iy_last, iy_last_p_o, MAX(iy_last, iy_last_p_o) + blk_y[scale_ind], left_over_y[scale_ind]);*/

	    	/*new_cols[scale_ind] = floor((right[scale_ind] - left[scale_ind]) / blk_y[scale_ind]);
	    	new_rows[scale_ind] = floor((bottom[scale_ind] - top[scale_ind]) / blk_x[scale_ind]);
		Nlines_x[scale_ind] = new_rows[scale_ind] + (p_o - 1) * (new_rows[scale_ind] - 1);
		Nlines_y[scale_ind] = new_cols[scale_ind] + (p_o - 1) * (new_cols[scale_ind] - 1);

		left_over_x[scale_ind] = 0;
		left_over_y[scale_ind] = 0;*/

		/*printf("new_rows: %i\n", new_rows[scale_ind]);
		printf("new_cols: %i\n", new_cols[scale_ind]);
		printf("Nlines_x: %i\n", Nlines_x[scale_ind]);
		printf("Nlines_y: %i\n", Nlines_y[scale_ind]);*/	

	    	Nblocks[scale_ind] = Nlines_x[scale_ind] * Nlines_y[scale_ind];
		w_now = (w_now + 1) / 2;
		h_now = (h_now + 1) / 2;

	};

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
		else {
			// compute
			if ((ret = compute_vif(ref_diff_buf, dis_diff_buf, w, h, stride, stride, &score, &score_num, &score_den, scores, frm_idx, Nscales, blk_x, blk_y, new_rows, new_cols, Nlines_x, Nlines_y, p_o, left_over_x, left_over_y, is_it_diff)))
			{
			    printf("error: compute_vif_diff failed.\n");
			    fflush(stdout);
			    goto fail_or_end;
			}
		}
	}

	else {
		// compute
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
		fflush(stdout);
		printf("vif_den_diff: %d %f\n", frm_idx, score_den);
		fflush(stdout);
		for(int scale=0;scale<4;scale++){
		    printf("vif_num_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale]);
		    printf("vif_den_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}

	else { 
		printf("vif_num: %d %f\n", frm_idx, score_num);
		fflush(stdout);
		printf("vif_den: %d %f\n", frm_idx, score_den);
		fflush(stdout);
		for(int scale=0;scale<Nscales;scale++){
		    printf("vif_num_scale%d: %d %f\n", scale, frm_idx, scores[2*scale]);
		    printf("vif_den_scale%d: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
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

    aligned_free(prev_ref_buf);
    aligned_free(prev_dis_buf);

    aligned_free(temp_buf);

    return ret;
}
