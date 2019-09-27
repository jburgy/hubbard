function multiply(n, k, nCk, pairs, t, U, x) {
    const y = [];
    x.forEach((w, i) => {
        if (!w) {
            return;
        }
        const h = (i / nCk) | 0; // MSB
        const l = i - h * nCk;   // LSB
        const u = mask(h, n, k, nCk);
        const d = mask(l, n, k, nCk);

        y[i] = (y[i] || 0) + U * w * count(u & d);

        const tw = t * w;
        pairs.forEach((p) => {
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