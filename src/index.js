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

const cache = {};
const steps = 100;

function parameters(document) {
    return Object.fromEntries(Array.from(document
        .querySelectorAll('fieldset#model > input'), ({ id, valueAsNumber }) => [id, valueAsNumber]));
}

function geometry(document) {
    const { valueAsNumber: k } = document.querySelector('fieldset#geometry > input#k');
    const { list, value, valueAsNumber: n } = document.querySelector('fieldset#geometry > input#n');
    const { dataset: { u, v } } = list.querySelector(`option[value="${value}"]`);

    return { u: parseInt(u), v: parseInt(v), n, k };
}

export default function callback(document) {
    const { n, k, u, v } = geometry(document);
    const key = (n << 8) | k;
    if (!(key in cache)) {
        cache[key] = hubbard(n, k, u, v);
    }
    return cache[key].then((multiply) => {
        const { t, U } = parameters(document);
        const nCk = binom(n, k);
        const dimension = nCk * nCk;
        const buffer = lanczos(dimension, steps, (x, y) => multiply(t, U, x, y))

        const d = new Float64Array(buffer, 0, steps);
        const e = new Float64Array(buffer, steps * Float64Array.BYTES_PER_ELEMENT, steps);
        tql1(d, e);
        return d;
    });
}