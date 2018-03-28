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
#include <string.h>

#include "common_christos/alloc.h"
#include "common_christos/tools.h"
#include "common_christos/file_io.h"
#include "adm_options.h"
#include "adm_tools.h"

#ifdef ADM_OPT_SINGLE_PRECISION
  typedef float number_t;
  typedef adm_dwt_band_t_s adm_dwt_band_t;

  #define read_image_b  read_image_b2s
  #define read_image_w  read_image_w2s
  #define adm_dwt2      adm_dwt2_s
  #define adm_decouple  adm_decouple_s
  #define adm_csf       adm_csf_s
  #define adm_cm_thresh adm_cm_thresh_s
  #define adm_cm        adm_cm_s
  #define adm_sum_cube  adm_sum_cube_s
  #define adm_sum_cube_dump_blocks  adm_sum_cube_s_dump_blocks
#else
  typedef double number_t;
  typedef adm_dwt_band_t_d adm_dwt_band_t;

  #define read_image_b  read_image_b2d
  #define read_image_w  read_image_w2d
  #define adm_dwt2      adm_dwt2_d
  #define adm_decouple  adm_decouple_d
  #define adm_csf       adm_csf_d
  #define adm_cm_thresh adm_cm_thresh_d
  #define adm_cm        adm_cm_d
  #define adm_sum_cube  adm_sum_cube_d
#endif

static char *init_dwt_band(adm_dwt_band_t *band, char *data_top, size_t buf_sz_one)
{
    band->band_a = (number_t *)data_top; data_top += buf_sz_one;
    band->band_h = (number_t *)data_top; data_top += buf_sz_one;
    band->band_v = (number_t *)data_top; data_top += buf_sz_one;
    band->band_d = (number_t *)data_top; data_top += buf_sz_one;
    return data_top;
}

int get_block(const number_t *ref, int ref_stride, int ref_stride_block, int start_x, int roi_x, int start_y, int roi_y, number_t ref_block[roi_x*roi_y])
{

    int ref_stride_ = ref_stride / sizeof(number_t);
    int ref_stride_block_ = ref_stride_block / sizeof(number_t);

    for (int i = start_x; i < start_x + roi_x; ++i)
    {
        for (int j = start_y; j < start_y + roi_y; ++j)
        {
	    ref_block[(i - start_x) * ref_stride_block_ + j - start_y] = ref[i * ref_stride_ + j];
	    /*printf("%f %f\n", ref_block[(i - start_x) * ref_stride_block_ + j - start_y], ref[i * ref_stride_ + j]);*/
	    /*printf("%f ", ref_block[(i - start_x) * ref_stride_block_ + j - start_y]);*/
        }
	/*printf("\n");*/
    }
    
    return 0;
}

int compute_adm(const number_t *ref, const number_t *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor, int Nscales, int blk_x[Nscales], int blk_y[Nscales], int Nblocks[Nscales], int px_stride[Nscales], int left[Nscales], int right[Nscales], int top[Nscales], int bottom[Nscales], int new_rows[Nscales], int new_cols[Nscales], int sum_of_Nblocks[Nscales], int frm_idx, int p_o, int Nlines_x[Nscales], int Nlines_y[Nscales], float adm_num_blk[Nscales], float adm_den_blk[Nscales], int left_over_x[Nscales], int left_over_y[Nscales], int is_it_diff)
{
#ifdef ADM_OPT_SINGLE_PRECISION
    double numden_limit = 1e-2 * (w * h) / (1920.0 * 1080.0);
#else
    double numden_limit = 1e-10 * (w * h) / (1920.0 * 1080.0);
#endif
    number_t *data_buf = 0;
    char *data_top;

    number_t *ref_scale;
    number_t *dis_scale;

    adm_dwt_band_t ref_dwt2;
    adm_dwt_band_t dis_dwt2;

    adm_dwt_band_t decouple_r;
    adm_dwt_band_t decouple_a;

    adm_dwt_band_t csf_o;
    adm_dwt_band_t csf_r;
    adm_dwt_band_t csf_a;

    number_t *mta;

    adm_dwt_band_t cm_r;

    const number_t *curr_ref_scale = ref;
    const number_t *curr_dis_scale = dis;
    int curr_ref_stride = ref_stride;
    int curr_dis_stride = dis_stride;

    int orig_h = h;

    int buf_stride = ALIGN_CEIL(((w + 1) / 2) * sizeof(number_t));
    size_t buf_sz_one = (size_t)buf_stride * ((h + 1) / 2);

    double num = 0;
    double den = 0;

    int scale;
    int ret = 1;

    if (SIZE_MAX / buf_sz_one < 35)
    {
        printf("error: SIZE_MAX / buf_sz_one < 35, buf_sz_one = %lu.\n", buf_sz_one);
        fflush(stdout);
        goto fail;
    }

    if (!(data_buf = aligned_malloc(buf_sz_one * 35, MAX_ALIGN)))
    {
        printf("error: aligned_malloc failed for data_buf.\n");
        fflush(stdout);
        goto fail;
    }

    data_top = (char *)data_buf;

    ref_scale = (number_t *)data_top; data_top += buf_sz_one;
    dis_scale = (number_t *)data_top; data_top += buf_sz_one;

    data_top = init_dwt_band(&ref_dwt2, data_top, buf_sz_one);
    data_top = init_dwt_band(&dis_dwt2, data_top, buf_sz_one);
    data_top = init_dwt_band(&decouple_r, data_top, buf_sz_one);
    data_top = init_dwt_band(&decouple_a, data_top, buf_sz_one);
    data_top = init_dwt_band(&csf_o, data_top, buf_sz_one);
    data_top = init_dwt_band(&csf_r, data_top, buf_sz_one);
    data_top = init_dwt_band(&csf_a, data_top, buf_sz_one);

    mta = (number_t *)data_top; data_top += buf_sz_one;

    data_top = init_dwt_band(&cm_r, data_top, buf_sz_one);

    float num_blk_sm = 0.0f;
    float den_blk_sm = 0.0f;

    /*printf("Nlinesx: %i\n", Nlines_x[scale]);
    float adm_num_blk[2] = 0.0f;*/
    /*float adm_den_blk[Nblocks0];
    for (int i = 0; i < Nblocks[0]; ++i) {
	adm_num_blk[i] = 0.0;
	adm_den_blk[i] = 0.0;
    }*/

    for (scale = 0; scale < Nscales; ++scale) {
#ifdef ADM_OPT_DEBUG_DUMP
        char pathbuf[256];
#endif
        float num_scale = 0.0;
        float den_scale = 0.0;
	float num_now, den_now;

        adm_dwt2(curr_ref_scale, &ref_dwt2, w, h, curr_ref_stride, buf_stride);
        adm_dwt2(curr_dis_scale, &dis_dwt2, w, h, curr_dis_stride, buf_stride);

        w = (w + 1) / 2;
        h = (h + 1) / 2;

        adm_decouple(&ref_dwt2, &dis_dwt2, &decouple_r, &decouple_a, w, h, buf_stride, buf_stride, buf_stride, buf_stride);

        adm_csf(&ref_dwt2, &csf_o, orig_h, scale, w, h, buf_stride, buf_stride);
        adm_csf(&decouple_r, &csf_r, orig_h, scale, w, h, buf_stride, buf_stride);
        adm_csf(&decouple_a, &csf_a, orig_h, scale, w, h, buf_stride, buf_stride);

        adm_cm_thresh(&csf_a, mta, w, h, buf_stride, buf_stride);
        adm_cm(&csf_r, &cm_r, mta, w, h, buf_stride, buf_stride, buf_stride);

#ifdef ADM_OPT_DEBUG_DUMP
        sprintf(pathbuf, "stage/ref[%d]_a.yuv", scale);
        write_image(pathbuf, ref_dwt2.band_a, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/ref[%d]_h.yuv", scale);
        write_image(pathbuf, ref_dwt2.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/ref[%d]_v.yuv", scale);
        write_image(pathbuf, ref_dwt2.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/ref[%d]_d.yuv", scale);
        write_image(pathbuf, ref_dwt2.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis[%d]_a.yuv", scale);
        write_image(pathbuf, dis_dwt2.band_a, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis[%d]_h.yuv", scale);
        write_image(pathbuf, dis_dwt2.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis[%d]_v.yuv", scale);
        write_image(pathbuf, dis_dwt2.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/dis[%d]_d.yuv", scale);
        write_image(pathbuf, dis_dwt2.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/r[%d]_h.yuv", scale);
        write_image(pathbuf, decouple_r.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/r[%d]_v.yuv", scale);
        write_image(pathbuf, decouple_r.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/r[%d]_d.yuv", scale);
        write_image(pathbuf, decouple_r.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/a[%d]_h.yuv", scale);
        write_image(pathbuf, decouple_a.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/a[%d]_v.yuv", scale);
        write_image(pathbuf, decouple_a.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/a[%d]_d.yuv", scale);
        write_image(pathbuf, decouple_a.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_o[%d]_h.yuv", scale);
        write_image(pathbuf, csf_o.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_o[%d]_v.yuv", scale);
        write_image(pathbuf, csf_o.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_o[%d]_d.yuv", scale);
        write_image(pathbuf, csf_o.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_r[%d]_h.yuv", scale);
        write_image(pathbuf, csf_r.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_r[%d]_v.yuv", scale);
        write_image(pathbuf, csf_r.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_r[%d]_d.yuv", scale);
        write_image(pathbuf, csf_r.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_a[%d]_h.yuv", scale);
        write_image(pathbuf, csf_a.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_a[%d]_v.yuv", scale);
        write_image(pathbuf, csf_a.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/csf_a[%d]_d.yuv", scale);
        write_image(pathbuf, csf_a.band_d, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/mta[%d].yuv", scale);
        write_image(pathbuf, mta, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/cm_r[%d]_h.yuv", scale);
        write_image(pathbuf, cm_r.band_h, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/cm_r[%d]_v.yuv", scale);
        write_image(pathbuf, cm_r.band_v, w, h, buf_stride, sizeof(number_t));

        sprintf(pathbuf, "stage/cm_r[%d]_d.yuv", scale);
        write_image(pathbuf, cm_r.band_d, w, h, buf_stride, sizeof(number_t));
#endif

	float num_scale_blocks_rowh[Nlines_x[scale] * Nlines_y[scale]];
	float den_scale_blocks_rowh[Nlines_x[scale] * Nlines_y[scale]];
	float num_scale_blocks_rowv[Nlines_x[scale] * Nlines_y[scale]];
	float den_scale_blocks_rowv[Nlines_x[scale] * Nlines_y[scale]];
	float num_scale_blocks_rowd[Nlines_x[scale] * Nlines_y[scale]];
	float den_scale_blocks_rowd[Nlines_x[scale] * Nlines_y[scale]];

	float num_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];
	float den_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];
	float adm_scale_blocks_row[Nlines_x[scale] * Nlines_y[scale]];

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) { 
		num_scale_blocks_rowh[i] = 0.0f; 
		num_scale_blocks_rowv[i] = 0.0f; 
		num_scale_blocks_rowd[i] = 0.0f;
		num_scale_blocks_row[i] = 0.0f;
		den_scale_blocks_rowh[i] = 0.0f; 
		den_scale_blocks_rowv[i] = 0.0f; 
		den_scale_blocks_rowd[i] = 0.0f;
		den_scale_blocks_row[i] = 0.0f;
		adm_scale_blocks_row[i] = 0.0f;
	}

	/*printf("%f\n", accum_blocks_row[0]); free(accum_blocks_row);
	printf("scale %i, Nblocks: %i\n", scale, Nblocks[scale]);
        printf("width: %i height %i rows: %i cols: %i right: %i left: %i top: %i bottom: %i \n", w, h, new_rows[scale], new_cols[scale], right[scale], left[scale], top[scale], bottom[scale]);*/
	
	/*if (scale == 3) {*/
	/*printf("Nlines_x out: %i\n", Nlines_x[scale]);*/
        adm_sum_cube_s_dump_blocks(cm_r.band_h, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], num_scale_blocks_rowh, left_over_x[scale], left_over_y[scale]);
        adm_sum_cube_s_dump_blocks(csf_o.band_h, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], den_scale_blocks_rowh, left_over_x[scale], left_over_y[scale]);
	/*}*/

        num_scale += adm_sum_cube(cm_r.band_h, w, h, buf_stride, border_factor);
        den_scale += adm_sum_cube(csf_o.band_h, w, h, buf_stride, border_factor);

	/*if (scale == 3) { 
		printf("Nblocks: %i\n", Nblocks[scale]);
		printf("num scale: %f\n", num_scale);
		printf("num scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", num_scale_blocks_row[i]); };
		printf("\nden scale: %f\n", den_scale);
		printf("den scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", den_scale_blocks_row[i]); };
		printf("\n");
	}

	printf("/////////////////////////////////////////\n"); */

	adm_sum_cube_s_dump_blocks(cm_r.band_v, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], num_scale_blocks_rowv, left_over_x[scale], left_over_y[scale]);
	adm_sum_cube_s_dump_blocks(csf_o.band_v, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], den_scale_blocks_rowv, left_over_x[scale], left_over_y[scale]);

        num_scale += adm_sum_cube(cm_r.band_v, w, h, buf_stride, border_factor);
        den_scale += adm_sum_cube(csf_o.band_v, w, h, buf_stride, border_factor);

	/*if (scale == 3) { 
		printf("Nblocks: %i\n", Nblocks[scale]);
		printf("num scale: %f\n", num_scale);
		printf("num scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", num_scale_blocks_row[i]); };
		printf("\nden scale: %f\n", den_scale);
		printf("den scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", den_scale_blocks_row[i]); };
		printf("\n");
	}*/

	/*printf("/////////////////////////////////////////\n");*/

	adm_sum_cube_s_dump_blocks(cm_r.band_d, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], num_scale_blocks_rowd, left_over_x[scale], left_over_y[scale]);
	adm_sum_cube_s_dump_blocks(csf_o.band_d, w, h, buf_stride, blk_x[scale], blk_y[scale], new_cols[scale], new_rows[scale], left[scale], right[scale], top[scale], bottom[scale], p_o, Nlines_x[scale], Nlines_y[scale], den_scale_blocks_rowd, left_over_x[scale], left_over_y[scale]);

        num_scale += adm_sum_cube(cm_r.band_d, w, h, buf_stride, border_factor);
        den_scale += adm_sum_cube(csf_o.band_d, w, h, buf_stride, border_factor);

        num += num_scale;
        den += den_scale;

	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) { 
		num_scale_blocks_row[i] = num_scale_blocks_rowh[i] + num_scale_blocks_rowv[i] + num_scale_blocks_rowd[i];
	}
	/*printf("Length of h: %i, Length of d:%p\n", sizeof(den_scale_blocks_rowh), den_scale_blocks_rowd);*/
	for (int i = 0; i < Nlines_x[scale] * Nlines_y[scale]; i++) { 
		den_scale_blocks_row[i] = den_scale_blocks_rowh[i] + den_scale_blocks_rowv[i] + den_scale_blocks_rowd[i]; 
	}

	/*if (scale == 3) { 
		printf("Nblocks: %i\n", Nblocks[scale]);
		printf("num scale: %f\n", num_scale);
		printf("num scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", num_scale_blocks_row[i]); };
		printf("\nden scale: %f\n", den_scale);
		printf("den scale blocks: ");
		for (int i = 0; i < Nblocks[scale]; i++) { printf("%f ", den_scale_blocks_row[i]); };
		printf("\n");
	}*/

        /* Copy DWT2 approximation band to buffer for next scale. */
        adm_buffer_copy(ref_dwt2.band_a, ref_scale, w * sizeof(number_t), h, buf_stride, buf_stride);
        adm_buffer_copy(dis_dwt2.band_a, dis_scale, w * sizeof(number_t), h, buf_stride, buf_stride);

        curr_ref_scale = ref_scale;
        curr_dis_scale = dis_scale;
        curr_ref_stride = buf_stride;
        curr_dis_stride = buf_stride;
#ifdef ADM_OPT_DEBUG_DUMP
        printf("num: %f\n", num);
        printf("den: %f\n", den);
#endif
        scores[2*scale+0] = num_scale;
        scores[2*scale+1] = den_scale;

	/***************************************************/

	/*if (scale > 0) {

		/*printf("for num: from %i to %i\n", sum_of_Nblocks[scale-1], sum_of_Nblocks[scale-1] + Nblocks[scale] - 1);
		printf("for den: from %i to %i\n", sum_of_Nblocks[scale-1] + Nblocks[scale], sum_of_Nblocks[scale-1] + 2*Nblocks[scale] - 1);
		for (int i = 0; i < Nblocks[scale]; i++)
			{
				scores_block[sum_of_Nblocks[scale-1] + i] = num_scale_blocks_row[i];
				scores_block[sum_of_Nblocks[scale-1] + Nblocks[scale] + i] = den_scale_blocks_row[i];
			}
	}
	else {

		/*printf("for num: from %i to %i\n", 0, Nblocks[scale] - 1);
		printf("for den: from %i to %i\n", Nblocks[scale], 2*Nblocks[scale] - 1);
		for (int i = 0; i < Nblocks[scale]; i++)
			{
				scores_block[i] = num_scale_blocks_row[i];
				scores_block[Nblocks[scale] + i] = den_scale_blocks_row[i];
			}
	}*/

	for (int i = 0; i < Nblocks[scale]; i++) {
		num_now = num_scale_blocks_row[i];
		num_now = num_now < numden_limit ? 0 : num_now;
		den_now = den_scale_blocks_row[i];
		den_now = den_now < numden_limit ? 0 : den_now;
		num_scale_blocks_row[i] = num_now;
		den_scale_blocks_row[i] = den_now;
		if (den_now == 0.0) {
			adm_scale_blocks_row[i] = 1.0f;
		}
		else {
			adm_scale_blocks_row[i] = num_scale_blocks_row[i] / den_scale_blocks_row[i];
		}
	}

	/*printf("adm, blk_x %i, blk_y %i, width %i, height %i\n", blk_x[scale], blk_y[scale], w, h);*/

	/* assumes that all scales have the same number of blocks */
	for (int i = 0; i < Nblocks[0]; i++) {
		adm_num_blk[i] += num_scale_blocks_row[i];
		adm_den_blk[i] += den_scale_blocks_row[i];
	}

	if (is_it_diff == 1) {
		printf("adm_num_scale%i_diff_blk: %i", scale, frm_idx); 
	}
	else {
		printf("adm_num_scale%i_blk: %i", scale, frm_idx); 
	}
	
	for (int i = 0; i < Nblocks[scale]; i++) {
		printf(" %f", num_scale_blocks_row[i]);
	}
	printf("\n");
        fflush(stdout);

	if (is_it_diff == 1) {
		printf("adm_den_scale%i_diff_blk: %i", scale, frm_idx); 
	}
	else {
		printf("adm_den_scale%i_blk: %i", scale, frm_idx); 
	}
	for (int i = 0; i < Nblocks[scale]; i++) {
		printf(" %f", den_scale_blocks_row[i]);
	}
	printf("\n");
        fflush(stdout);

	if (is_it_diff == 1) {
		printf("adm_scale%i_diff_blk: %i", scale, frm_idx); 
	}
	else {
		printf("adm_scale%i_blk: %i", scale, frm_idx); 
	}
	for (int i = 0; i < Nblocks[scale]; i++) {
		printf(" %f", adm_scale_blocks_row[i]);
	}
	printf("\n");
        fflush(stdout);

	for (int i = 0; i < Nblocks[scale]; i++) {

		num_blk_sm += num_scale_blocks_row[i];
		den_blk_sm += den_scale_blocks_row[i];

	}

    }

    if (is_it_diff == 1) {
    	printf("adm_num_diff_blk: %i", frm_idx);
    }
    else {
    	printf("adm_num_blk: %i", frm_idx);
    }
    for (int i = 0; i < Nblocks[0]; i++) {
    	printf(" %f", adm_num_blk[i]);
    }
    printf("\n");
    fflush(stdout);

    if (is_it_diff == 1) {
    	printf("adm_den_diff_blk: %i", frm_idx);
    }
    else {
    	printf("adm_den_blk: %i", frm_idx);
    }
    for (int i = 0; i < Nblocks[0]; i++) {
    	printf(" %f", adm_den_blk[i]);
    }
    printf("\n");
    fflush(stdout);
    
    if (is_it_diff == 1) {
    	printf("adm_sm_diff_blk: %i %f\n", frm_idx, num_blk_sm / den_blk_sm); 
    }
    else {
	printf("adm_sm_blk: %i %f\n", frm_idx, num_blk_sm / den_blk_sm); 
    }
    fflush(stdout);

  /*    for (int i = 0; i < scores_block_length/2; i++) {
	/* can belong to a scale or a block within a scale
	
	if (den_now == 0.0) {
		scores_block[i] = 1.0d;
	}
	else {
		scores_block[i] = num_now / den_now;
	}
	scores_block[i] = num_now;
	scores_block[i] = den_now;

    }*/

    num = num < numden_limit ? 0 : num;
    den = den < numden_limit ? 0 : den;

    if (den == 0.0)
    {
        *score = 1.0f;
    }
    else
    {
        *score = num / den;
    }
    *score_num = num;
    *score_den = den;

    ret = 0;

fail:
    aligned_free(data_buf);
    return ret;
}

int adm(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int start_x, int start_y, int roi_x, int roi_y, int blk_x_val, int blk_y_val, int p_o, int is_it_diff)
{
    double score = 0;
    double score_num = 0;
    double score_den = 0;
    double scores[2*4];

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

    int Nscales = 4;

    int px_stride[Nscales], left[Nscales], right[Nscales], top[Nscales], bottom[Nscales], Nblocks[Nscales];	
    int new_rows[Nscales], new_cols[Nscales], blk_x[Nscales], blk_y[Nscales], sum_of_Nblocks[Nscales];
    int Nlines_x[Nscales], Nlines_y[Nscales], left_over_x[Nscales], left_over_y[Nscales];
    float num_now, den_now;
    int i, ix_last, iy_last, ix_last_p_o, iy_last_p_o;
    int w_now, h_now;

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
        printf("error: fopen ref_path %s failed\n", ref_path);
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

        // compute

	int stride_block = ALIGN_CEIL(roi_y * sizeof(number_t));
	w_now = w;
        h_now = h;
	number_t ref_block[roi_x*roi_y];
	number_t dis_block[roi_x*roi_y];

	/*printf("Reference block #%d: \n", frm_idx);*/
	/*get_block(ref_buf, stride, stride_block, start_x, roi_x, start_y, roi_y, ref_block);*/
	
	/*printf("Distorted block #%d: \n", frm_idx);*/
	/*get_block(dis_buf, stride, stride_block, start_x, roi_x, start_y, roi_y, dis_block);*/

	int buf_stride = ALIGN_CEIL(((w + 1) / 2) * sizeof(number_t));

	for (int scale_ind = 0; scale_ind < Nscales; scale_ind++) { sum_of_Nblocks[scale_ind] = 0; };

	for (int scale_ind = 0; scale_ind < Nscales; scale_ind++) { 
		
		w_now = (w_now + 1) / 2;
		h_now = (h_now + 1) / 2;
		Nblocks[scale_ind] = 0; 
		blk_x[scale_ind] = 0;
		blk_y[scale_ind] = 0;

		px_stride[scale_ind] = buf_stride / sizeof(number_t);
		/*left[scale_ind]   = w_now * ADM_BORDER_FACTOR - 0.5;
	    	top[scale_ind]    = h_now * ADM_BORDER_FACTOR - 0.5;*/
		left[scale_ind]   = MAX(w_now * ADM_BORDER_FACTOR - 0.5, 0); 
	    	top[scale_ind]    = MAX(h_now * ADM_BORDER_FACTOR - 0.5, 0);
		right[scale_ind]  = w_now - left[scale_ind];
		bottom[scale_ind] = h_now - top[scale_ind];  
		if ((blk_x_val == -1) || (blk_x_val == -1)) {
			blk_x[scale_ind] = bottom[scale_ind] - top[scale_ind];
			blk_y[scale_ind] = right[scale_ind] - left[scale_ind];
    		}
		else {
			if (scale_ind == 0) {
				/* remember that adm starts of at the lower scale */
				blk_x[scale_ind] = MIN((blk_x_val + 1) / 2, bottom[scale_ind] - top[scale_ind]);
				blk_y[scale_ind] = MIN((blk_y_val + 1) / 2, right[scale_ind] - left[scale_ind]);
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
		Nlines_y[scale_ind] = new_cols[scale_ind] + (p_o - 1) * (new_cols[scale_ind] - 1);*/

	    	Nblocks[scale_ind] = Nlines_x[scale_ind] * Nlines_y[scale_ind];

		/*printf("new_rows: %i", new_rows[scale_ind]);
		printf(" new_cols: %i", new_cols[scale_ind]);
		printf(" Nlines_x: %i", Nlines_x[scale_ind]);
		printf(" Nlines_y: %i\n", Nlines_y[scale_ind]);*/

		/*if (scale_ind == 0) {
			sum_of_Nblocks[scale_ind] = 2*Nblocks[scale_ind];
		}
		else {
			for (int scale_ind_inner = 0; scale_ind_inner < scale_ind; scale_ind_inner++) {
				sum_of_Nblocks[scale_ind] += sum_of_Nblocks[scale_ind_inner];
			};
			sum_of_Nblocks[scale_ind] += 2*Nblocks[scale_ind];
		}*/

	};

	float adm_num_blk[Nblocks[0]], adm_den_blk[Nblocks[0]];
    	for (int i = 0; i < Nblocks[0]; ++i) {
		adm_num_blk[i] = 0.0;
		adm_den_blk[i] = 0.0;
    	}

	/*for (int i = 0; i < Nscales; i++) { printf("Nblocks %i Nlines_x %i Nlines_y %i new_cols %i new_rows %i ", Nblocks[i], Nlines_x[i], Nlines_y[i], new_cols[i], new_rows[i]); } ;
	printf("\n");*/
	/*printf("Nblocks")*/
	/*printf("width %i", w);*/

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
			ret = compute_adm(ref_diff_buf, dis_diff_buf, w, h, stride_block, stride_block, &score, &score_num, &score_den, scores, ADM_BORDER_FACTOR, Nscales, blk_x, blk_y, Nblocks, px_stride, left, right, top, bottom, new_rows, new_cols, sum_of_Nblocks, frm_idx, p_o, Nlines_x, Nlines_y, adm_num_blk, adm_den_blk, left_over_x, left_over_y, is_it_diff);
			/*ret = compute_adm(ref_buf, dis_buf, w, h, stride, stride, &score, &score_num, &score_den, scores, ADM_BORDER_FACTOR)*/
			if ((ret))
			{
			    printf("error: compute_adm_diff failed.\n");
			    fflush(stdout);
			    goto fail_or_end;
			}
		}
	}

	else {
		ret = compute_adm(ref_buf, dis_buf, w, h, stride_block, stride_block, &score, &score_num, &score_den, scores, ADM_BORDER_FACTOR, Nscales, blk_x, blk_y, Nblocks, px_stride, left, right, top, bottom, new_rows, new_cols, sum_of_Nblocks, frm_idx, p_o, Nlines_x, Nlines_y, adm_num_blk, adm_den_blk, left_over_x, left_over_y, is_it_diff);
		if ((ret))
		{
		    printf("error: compute_adm failed.\n");
		    fflush(stdout);
		    goto fail_or_end;
		}
	}

        // print

	if (is_it_diff == 1) {
		printf("adm_diff: %d %f\n", frm_idx, score);
		fflush(stdout);
		printf("adm_num_diff: %d %f\n", frm_idx, score_num);
		fflush(stdout);
		printf("adm_den_diff: %d %f\n", frm_idx, score_den);
		fflush(stdout);
		for(int scale=0;scale<4;scale++){
		    printf("adm_num_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale]);
		    printf("adm_den_scale%d_diff: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}

	else {
		printf("adm: %d %f\n", frm_idx, score);
		fflush(stdout);
		printf("adm_num: %d %f\n", frm_idx, score_num);
		fflush(stdout);
		printf("adm_den: %d %f\n", frm_idx, score_den);
		fflush(stdout);
		for(int scale=0;scale<4;scale++){
		    printf("adm_num_scale%d: %d %f\n", scale, frm_idx, scores[2*scale]);
		    printf("adm_den_scale%d: %d %f\n", scale, frm_idx, scores[2*scale+1]);
		}
	}

        /* remember to print for frame 0 as well */
        if (is_it_diff == 1 && frm_idx == 0) {

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
