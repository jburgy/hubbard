import hubbard from './multiply.js';
import lanczos from './lanczos.js';
import tql1 from './tql1.js';

function binom(n, k) {
    let res = 1;
    for (let t = Math.min(k, n - k); t >= 1; --t) {
        res *= n-- / t;
    }
    return res | 0;
}

export default function diagonalize(N, K, u, v, t, U, steps) {
    return hubbard(N, K, u, v)
        .then((multiply) => {
            const nCk = binom(N, K);
            const dimension = nCk * nCk;
            const buffer = lanczos(dimension, steps, (x, y) => multiply(t, U, x, y))

            const d = new Float64Array(buffer, 0, steps);
            const e = new Float64Array(buffer, steps * Float64Array.BYTES_PER_ELEMENT, steps);
            tql1(d, e);
            return d;
        });
}