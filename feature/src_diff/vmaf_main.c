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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common_christos/cpu.h"

int adm(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int start_x, int start_y, int roi_x, int roi_y, int blk_x_val, int blk_y_val, int overlap, int is_it_diff);
int ansnr(const char *ref_path, const char *dis_path, int w, int h, const char *fmt);
int vif(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int blk_x_val, int blk_y_val, int overlap, int is_it_diff);
int motion(const char *dis_path, int w, int h, const char *fmt, int blk_x_val, int blk_y_val, int overlap);
int all(const char *ref_path, const char *dis_path, int w, int h, const char *fmt, int start_x, int start_y, int roi_x, int roi_y, int blk_x_val, int blk_y_val, int overlap, int is_it_diff);

enum vmaf_cpu cpu; // global

static void usage(void)
{
    puts("usage: vmaf app fmt ref dis w h [start_x start_y roi_x roi_y blk_x blk_y overlap is_it_diff]\n"
         "apps:\n"
         "\tadm\n"
         "\tansnr\n"
         "\tmotion\n"
         "\tvif\n"
         "\tall\n"
         "fmts:\n"
         "\tyuv420p\n"
         "\tyuv422p\n"
         "\tyuv444p\n"
         "\tyuv420p10le\n"
         "\tyuv422p10le\n"
         "\tyuv444p10le\n\n"
         "start_x: [0, h], default = 0\n"
         "start_y: [0, w], default = 0\n"
         "roi_x: [1, h], default = h\n"
         "roi_y: [1, w], default = w\n"
         "blk_x: 2^k for some k, default = -1 (whole frame)\n"
         "blk_y: 2^l for some l, default = -1 (whole frame)\n"
	 "overlap: some integer that makes sense e.g. 1 (no overlap) or 2 (50% per axis), default = 1\n"
	 "is_it_diff: 0 (applied on regular frames) or 1 (use frame diff), default = 0\n"
    );
}

int main(int argc, const char **argv)
{
    const char *app;
    const char *ref_path;
    const char *dis_path;
    const char *fmt;
    int w;
    int h;
    int ret;
    int start_x, start_y;
    int roi_x, roi_y;
    int blk_x = -1, blk_y = -1;
    int overlap = 1;
    int is_it_diff = 0;

    if (argc < 7) {
        usage();
        return 2;
    }

    app      = argv[1];
    fmt      = argv[2];
    ref_path = argv[3];
    dis_path = argv[4];
    w        = atoi(argv[5]);
    h        = atoi(argv[6]);

    if (w <= 0 || h <= 0) {
        usage();
        return 2;
    }

    if (argc > 10) {

        start_x  = atoi(argv[7]);
        start_y  = atoi(argv[8]);
        roi_x    = atoi(argv[9]);
        roi_y    = atoi(argv[10]);

	if (argc > 11) {
		blk_x = atoi(argv[11]);
		blk_y = atoi(argv[12]);
	}

	if (argc >= 14) {
		overlap = atoi(argv[13]);
	}

	if (argc >= 15) {
		is_it_diff = atoi(argv[14]);
	}

	if ((start_x + roi_x > h) || (start_y + roi_y > w)) {
		usage();
		return 2;
	}

    }

    else {

        start_x  = 0;
        start_y  = 0;
        roi_x    = h;
        roi_y    = w;

    }

    cpu = cpu_autodetect();

    if (!strcmp(app, "adm"))
        ret = adm(ref_path, dis_path, w, h, fmt, start_x, start_y, roi_x, roi_y, blk_x, blk_y, overlap, is_it_diff);
    else if (!strcmp(app, "ansnr"))
        ret = ansnr(ref_path, dis_path, w, h, fmt);
    else if (!strcmp(app, "vif"))
        ret = vif(ref_path, dis_path, w, h, fmt, blk_x, blk_y, overlap, is_it_diff);
    else if (!strcmp(app, "motion"))
        ret = motion(ref_path, w, h, fmt, blk_x, blk_y, overlap);
    else if (!strcmp(app, "all"))
        ret = all(ref_path, dis_path, w, h, fmt, start_x, start_y, roi_x, roi_y, blk_x, blk_y, overlap, is_it_diff);
    else
        return 2;

    if (ret)
        return ret;

    return 0;
}
