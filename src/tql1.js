/**
 * this subroutine is a translation of the algol procedure tql1,
 * num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
 * wilkinson.
 * handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
 * https://www.jstor.org/stable/2005222
 * 
 * this subroutine finds the eigenvalues of a symmetric
 * tridiagonal matrix by the ql method.
 * 
 * on input
 *  
 *    d contains the diagonal elements of the input matrix.
 *  
 *    e contains the super-diagonal elements of the input matrix.
 *  
 *  on output
 * 
 *    d contains the eigenvalues in ascending order.  if an
 *      error exit is made, the eigenvalues are correct and
 *      ordered for indices 0, 1, ...rc - 1, but may not be
 *      the smallest eigenvalues.
 * 
 *    e has been destroyed.
 * 
 *    result is set to
 *      zero       for normal return,
 *      j          if the j-th eigenvalue has not been
 *                  determined after 30 iterations.
 * 
 * questions and comments should be directed to burton s. garbow,
 * mathematics and computer science div, argonne national laboratory
 * 
 * this version dated august 1983.
 * 
 * @param {FloatArray64} d - diagonal elements of the input matrix
 * @param {FloatArray64} e - super-diagonal elements of the input matrix
 * @returns {number} m - number of eigenvalues which could be determined
 *                       in 30 iterations or fewer
 */
function tql1(d, e) {
    const { length: n } = d;
    let f = 0;
    let test = 0;

    e[n - 1] = 0;
    for (let l = 0; l < n; ++l) {
        const h = Math.abs(d[l]) + Math.abs(e[l]);
        if (test < h) {
            test = h;
        }
        // .......... look for small sub-diagonal element ..........
        let m = l;
        while (test + Math.abs(e[m]) !== test) {
            ++m;
            /* .......... e(n) is always zero, so there is no exit
                          through the bottom of the loop .......... */
        }
        let converged = m === l;
        for (let j = 0; !converged && j < 30; ++j) {
            // .......... form shift ..........
            const l1 = l + 1;
            const l2 = l1 + 1;
            let g = d[l]
            let p = (d[l1] - g) / (2 * e[l]);
            let r = p + Math.hypot(p, 1) * (p < 0 ? -1 : 1);
            d[l] = e[l] / r;
            d[l1] = e[l] * r;
            const dl1 = d[l1];
            let h = g - d[l];

            for (let i = l2; i < n; ++i) {
                d[i] -= h;
            }
            f += h;
            // .......... ql transformation ..........
            p = d[m];
            let c = 1;
            let c2 = 1;
            const el1 = e[l1];
            let s = 0;
            const mml = m - l

            let s2;
            let c3;
            for (let i = m - 1; i >= l; --i) {
                c3 = c2;
                c2 = c;
                s2 = s;
                g = c * e[i];
                h = c * p
                r = Math.hypot(p, e[i]);
                e[i + 1] = s * r;
                s = e[i] / r;
                c = p / r;
                p = c * d[i] - s * g;
                d[i + 1] = h + s * (c * g + s * d[i]);
            }

            p = -s * s2 * c3 * el1 * e[l] / dl1;
            e[l] = s * p;
            d[l] = c * p;
            converged = test + Math.abs(e[l]) === test;
        }
        if (!converged) {
            return l + 1;
        }
        
        const p = d[l] + f;
        // .......... order eigenvalues ..........
        let i = l;
        for (; i; --i) {
            if (p >= d[i - 1]) {
                break;
            }
            d[i] = d[i - 1];
        }
        d[i] = p;
    }
    return 0;
}
