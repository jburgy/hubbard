import neighbors from './lattice.js';

function multiply(nCk, mask, index, count, hops, t, U, x) {
    const y = [];
    x.forEach((w, i) => {
        if (!w) {
            return;
        }
        const h = (i / nCk) | 0; // MSB
        const l = i - h * nCk;   // LSB
        const u = mask(h);
        const d = mask(l);

        y[i] = (y[i] || 0) + U * w * count(u & d);

        const tw = t * w;
        hops.forEach((p) => {
            // fermions transfer when only one site is occupied
            const a = u & p;
            if (a && a !== p) {
                const j = index(u ^ p) * nCk + l;
                y[j] = (y[j] || 0) - tw;
            }
            const b = d & p;
            if (b && b !== p) {
                const j = h * nCk + index(d ^ p);
                y[j] = (y[j] || 0) - tw;
            }
        });
    });
    return y;
}

export default function hubbard(N, K, u, v) {
    return WebAssembly
        .instantiateStreaming(fetch('../out/combinations.wasm'))
        .then(({ instance: { exports: { index, mask, count } } }) => {
            const nCk = index(((1 << K) - 1) << (N - K)) + 1;
            const n = N - 1;
            const k = K;
            const mask$ = i => mask(i, n, k, nCk);

            const hops = neighbors(u, v)
                .flatMap((x, i) => x.filter((_, i) => !(i & 1)).map(j => (1 << i) | (1 << j)));

            return (t, U, x) => multiply(nCk, mask$, index, count, hops, t, U, x);
        });
}