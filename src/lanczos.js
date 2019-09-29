/**
 * Inner product
 * 
 * @param {FloatArray64} v 
 * @param {Float32Array} w 
 * @returns number
 */
function inner(v, w) {
    const { length } = v;
    let result = 0;

    for (let i = 0; i < length; ++i) {
        result += v[i] * w[i];
    }
    return result;
}

/**
 * @callback multiplyCallback
 * @param {Float64Array} x - input
 * @param {Float64Array} y - output
 */

/**
 * 
 * @param {number} dimension - size of Hermitian matrix
 * @param {number} steps - number of iterations
 * @param {multiplyCallback} multiply - v => A v
 * @returns {ArrayBuffer} diagonal and sub-diagonal concatenated into contiguous buffer
 */
export default function lanczos(dimension, steps, multiply) {
    // allocate all memory upfront for clarity
    const buffer = new ArrayBuffer(3 * dimension * Float64Array.BYTES_PER_ELEMENT);
    const result = new ArrayBuffer(2 * steps * Float64Array.BYTES_PER_ELEMENT);
    const alpha = new Float64Array(result, 0, steps);
    const beta = new Float64Array(result, steps * Float64Array.BYTES_PER_ELEMENT, steps);

    const offset = dimension * Float64Array.BYTES_PER_ELEMENT;
    const u = new Float64Array(buffer, 0 * offset, dimension); // v_{j - 1}
    const v = new Float64Array(buffer, 1 * offset, dimension); // v_j
    const w = new Float64Array(buffer, 2 * offset, dimension); // w_j
    
    v[0] = 1;
    multiply(v, w);

    const a = alpha[0] = inner(v, w);
    for (let i = 0; i < dimension; ++i) {
        w[i] -= a * v[i];
    }

    for (let step = 1; step < steps; ++step) {
        const b = beta[step - 1] = Math.sqrt(inner(w, w));

        for (let i = 0; i < dimension; ++i) {
            u[i] = v[i]; // save v_{j - 1} before overwriting
            v[i] = w[i] / b;
            w[i] = 0;
        }

        multiply(v, w);

        const a = alpha[step] = inner(v, w);
        for (let i = 0; i < dimension; ++i) {
            w[i] -= a * v[i] - b * u[i];
        }
    }

    return result;
}