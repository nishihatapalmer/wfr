/*
 * SMART: string matching algorithms research tool.
 * Copyright (C) 2012  Simone Faro and Thierry Lecroq
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 * contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
 * download the tool at: http://www.dmi.unict.it/~faro/smart/
 *
 *  THIS IS AN IMPLEMENTATION OF:
 *	ALGORITHM Weak Factor Recognizer (WFR) Using Q-Grams
 *  appeared in: Simone Faro, Domenico Cantone and Arianna Pavone.
 *  Speeding Up String Matching by Weak Factor Recognition.
 *  Proceedings of the Pague Stringology Conference 2017: pp.42-50
 *
 * PREPROCESSING:
 *		an hash value is computed for al factors of the pattern with length in [1..16]
 *		the computed hash value is always a number in [0...256*256]
 *		if w is a factor of x, and hash(w) is its hash value, than F[hash(w)]=TRUE, otherwise F[hash(w)]=FALSE
 * SEARCHING
 *		The algorithm searches for factors of the pattern using a weak recognition method
 *		the search phase is very similar to BOM.
 *		The window is scanned right to left and for each character a new hash value of the suffix of the window is computed.
 *		Let w be the suffix we scanned. If F[hash(w)]=TRUE we continue scanning the next character of the window.
 *		Otherwise we stop scanning (w is not a factor of the pattern) and jump to the right, like in BOM.
 */

#include "include/define.h"
#include "include/main.h"
#include "include/GRAPH.h"
#define Q 7
#define ALPHA 12

#define HASH(y, j) (((y)[(j)] << 12) + ((y)[(j) - 1] << 10) + ((y)[(j) - 2] << 8) + ((y)[(j) - 3] << 6) + ((y)[(j) - 4] << 4) + ((y)[(j) - 5] << 2) + (y)[(j) - 6])
#define ASIZE (1 << (ALPHA))
#define MASK ((ASIZE) - 1)

#define TEST_BIT(F, v) ((F)[((v) >> 5) & MASK] &  (1 << ((v) & 0x1F)))
#define SET_BIT(F, v)  ((F)[((v) >> 5) & MASK] |= (1 << ((v) & 0x1F)))

int preprocessing(unsigned char *x, int m, int *F) {
    int i, j, v;
    for (v = 0; v < ASIZE; v++) F[v] = 0;

    const int fact = m < 16 ? m : 16;
    for (i = m - 1; i >= Q - 1; i--)
    {
        const int stop = MAX(0, i - fact + 1);
        v = 0;
        for (j = i; j >= stop; j--) {
            v = (v << 2) + x[j];
            SET_BIT(F, v);
        }
    }
}

int search(unsigned char *x, int m, unsigned char *y, int n) {
    int i, j, k, count, test, h, F[ASIZE];
    if (m < Q) return -1;

    BEGIN_PREPROCESSING
    /* Preprocessing */
    const int plen = m;
    if (m % Q != 0) m -= (m % Q);
    const int mq = m - Q + 1;
    preprocessing(x, m, F);
    END_PREPROCESSING

    BEGIN_SEARCHING
    /* Searching */
    count = 0;
    j = m - 1;
    while (j < n) {

        h = HASH(y, j);
        i = j - m + Q;
        while((test = TEST_BIT(F, h)) && j > i + Q - 1) {
            j -= Q;
            h = (h << 14) + HASH(y, j);
        }

        if (j == i && test) {
            k = 0;
            i -= Q - 1;
            while (k < plen && x[k] == y[i + k]) k++;
            if (k == plen && i <= n - plen) count++;
        }

        j += mq;
    }
    END_SEARCHING
    return count;
}

