fetch('../out/main.wasm')
  .then(response => response.arrayBuffer())
  .then(bytes => WebAssembly.instantiate(bytes))
  .then(({ instance: { exports: { index, mask } } }) => {
    function generate(N, K) {
      const NCK = index(((1 << K) - 1) << (N - K)) + 1; // '1'.repeat(K) + '0'.repeat(N-K)
      const n = N - 1;
      const k = K;
      const nCk = NCK * (N - K) / N;
      const prefix = '0'.repeat(N - K);

      const list = [];
      for (let i = 0; i < NCK; ++i) {
        const m = mask(i, n, k, nCk);
        list.push(`${prefix}${m.toString(2)}`.slice(-N));
      }
      return list.join('<br>')
    }
    document.getElementById("container").innerHTML = generate(8, 5);
  })
  .catch(console.error);
