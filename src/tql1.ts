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
 * @param   {i32} d - diagonal elements of the input matrix
 * @param   {i32} e - super-diagonal elements of the input matrix
 * @returns {i32} m - number of eigenvalues which could be determined
 *                    in 30 iterations or fewer
 * 
 * d and e are pointers to the ```WebAssembly.Memory``` instance which
 * ```WebAssembly.instantiateStreaming``` passes in the .env.memory
 * component of its ```importObject``` second optional argument.
 * 
 * This requires using the ```--importMemory``` flag on the ```asc```
 * command line. Another recommended flag is ```--use Math=JSMath```.
 * 
 * Sample code:
 * :::javascript
 * const memory = new WebAssembly.Memory({ initial: 10 });
 *
 * WebAssembly.instantiateStreaming(fetch("tql1.wasm"), {
 *    env: { memory },
 *    Math,
 * })
 * .then(({ instance: { exports: { tql1 } } }) => {
 *    const size = Float64Array.BYTES_PER_ELEMENT;
 *    const n = 5;
 *    const d = new Float64Array(memory.buffer, 0, n);
 *    const e = new Float64Array(memory.buffer, n * size, n);
 *    for (let i = 0; i < n; ++i) {
 *        d[i] = 0;
 *        e[i] = 1;
 *    }
 *    console.assert(!tql1(0, n * size, n));
 *    console.log(d);
 * })
 * .catch(console.error);
 */
export function tql1(d: i32, e: i32, n: i32): i32 {
    let f:f64 = 0;
    let test:f64 = 0;
    const size = Float64Array.BYTES_PER_ELEMENT;
    const end = n * size;

    store<f64>(e + end - size, 0);
    for (let l = 0; l < end; l += size) {
        const h = Math.abs(load<f64>(d + l)) + Math.abs(load<f64>(e + l));
        if (test < h) {
            test = h;
        }
        // .......... look for small sub-diagonal element ..........
        let m = l;
        while (test + Math.abs(load<f64>(e + m)) !== test) {
            m += size;
            /* .......... e(n) is always zero, so there is no exit
                          through the bottom of the loop .......... */
        }
        let converged = m === l;
        for (let j = 0; !converged && j < 30; ++j) {
            // .......... form shift ..........
            const l1 = l + size;
            const l2 = l1 + size;
            let g = load<f64>(d + l);
            const el = load<f64>(e + l);
            let p = (load<f64>(d + l1) - g) / (2 * el);
            let r = p + Math.hypot(p, 1) * (p < 0 ? -1 : 1);
            store<f64>(d + l, el / r);
            store<f64>(d + l1, el * r);
            const dl1 = load<f64>(d + l1);
            let h = g - load<f64>(d + l);

            for (let i = l2; i < end; i += size) {
                store<f64>(d + i, load<f64>(d + i) - h);
            }
            f += h;
            // .......... ql transformation ..........
            p = load<f64>(d + m);
            let c:f64 = 1;
            let c2:f64 = 1;
            const el1 = load<f64>(e + l1);
            let s:f64 = 0;

            let s2:f64;
            let c3:f64;
            for (let i = m - size; i >= l; i -= size) {
                c3 = c2;
                c2 = c;
                s2 = s;
                const ei = load<f64>(e + i);
                g = c * ei;
                h = c * p
                r = Math.hypot(p, ei);
                store<f64>(e + i + size, s * r);
                s = ei / r;
                c = p / r;
                const di = load<f64>(d + i);
                p = c * di - s * g;
                store<f64>(d + i + size, h + s * (c * g + s * di));
            }

            p = -s * s2 * c3 * el1 * load<f64>(e + l) / dl1;
            store<f64>(e + l, s * p);
            store<f64>(d + l, c * p);
            converged = test + Math.abs(load<f64>(e + l)) === test;
        }
        if (!converged) {
            return l + 1;
        }
        
        const p = load<f64>(d + l) + f;
        // .......... order eigenvalues ..........
        let i = l;
        for (; i; i -= size) {
            const dim1 = load<f64>(d + i - size);
            if (p >= dim1) {
                break;
            }
            store<f64>(d + i, dim1);
        }
        store<f64>(d + i, p);
    }
    return 0;
}
