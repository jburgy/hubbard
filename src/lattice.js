const squares = Array.from({ length: 8 }, (_, i) => i * i);
export const sizes = Array
    .from({ length: 32 }, (_, i) => i + 1)
    .filter(n => squares.some(w => squares.includes(n - w)));

/**
 * Indices of 4 nearest neighbors of lattice contained in square
 * with base (u, v).
 * 
 * Coordinates (x, y) of points in square with base (u, v) satisfy
 *   0 <= u * x + v * y < u * u + v * v
 *   0 <= u * y - v * x < u * u + v * v
 * 
 * @param {number} n - number of sites
 * @returns {number[][]} -  result[i][j] is jth neighbor of ith site
 */
export default function neighbors(n) {
    const v = squares.findIndex(w => squares.includes(n - w));
    const u = squares.indexOf(n - v * v);
    const e = u * u + v * v;
    const reflect = (x, y) => {
        while (u * x + v * y < 0) {
            x += u;
            y += v;
        }
        while (u * x + v * y >= e) {
            x -= u;
            y -= v;
        }
        while (u * y - v * x < 0) {
            x -= v;
            y += u;
        }
        while (u * y - v * x >= e) {
            x += v;
            y -= u;
        }
        return [x, y]
    };
    const sites = [];
    const result = [];

    for (let y = 0; y < u + v; ++y) {
        const w = Math.ceil(y < u ?  -v * y / u : (u * y - e + 1) / v);
        for (let x = w; u * x + v * y < e && u * y - v * x >= 0; ++x) {
            sites.push([x, y]);
            result.push([
                reflect(x + 1, y),
                reflect(x - 1, y),
                reflect(x, y + 1),
                reflect(x, y - 1),
            ]);
        }
    }

    return result.map(a => a.map(i => sites.findIndex(j => !indexedDB.cmp(i, j))));
}
