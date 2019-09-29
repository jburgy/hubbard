before('Load WASM asynchronously', async function () {
    const { instance: { exports: { index, mask, count } } } = await WebAssembly
        .instantiateStreaming(fetch("../out/combinations.wasm"));

    this.count = count;
    this.generate = (N, K) => {
        const NCK = index(((1 << K) - 1) << (N - K)) + 1;
        const n = N - 1;
        const k = K;
        const nCk = NCK * (N - K) / N;

        return Array.from({ length: NCK }, (_, i) => mask(i, n, k, nCk));
    };
})

describe('Matches expectation', function () {
    it('Works for N=8, K=5', function () {
        const expected = [
            0b00011111, 0b00101111, 0b00110111, 0b00111011, 0b00111101, 0b00111110,
            0b01001111, 0b01010111, 0b01011011, 0b01011101, 0b01011110, 0b01100111,
            0b01101011, 0b01101101, 0b01101110, 0b01110011, 0b01110101, 0b01110110,
            0b01111001, 0b01111010, 0b01111100, 0b10001111, 0b10010111, 0b10011011,
            0b10011101, 0b10011110, 0b10100111, 0b10101011, 0b10101101, 0b10101110,
            0b10110011, 0b10110101, 0b10110110, 0b10111001, 0b10111010, 0b10111100,
            0b11000111, 0b11001011, 0b11001101, 0b11001110, 0b11010011, 0b11010101,
            0b11010110, 0b11011001, 0b11011010, 0b11011100, 0b11100011, 0b11100101,
            0b11100110, 0b11101001, 0b11101010, 0b11101100, 0b11110001, 0b11110010,
            0b11110100, 0b11111000,
        ];
        const realized = this.generate(8, 5);
        chai.expect(realized).to.deep.equal(expected);
    });

    it('count(0b10101010) === 4', function () {
        chai.expect(this.count(0b10101010)).to.be.equal(4);
    })
});