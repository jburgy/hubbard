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
        let z = u * x + v * y;
        while (z < 0) {
            x += u;
            y += v;
            z += e;
        }
        while (z >= e) {
            x -= u;
            y -= v;
            z -= e;
        }
        z = u * y - v * x;
        while (z < 0) {
            x -= v;
            y += u;
            z += e;
        }
        while (z >= e) {
            x += v;
            y -= u;
            z -= e;
        }
        return [x, y]
    };
    const sites = [];
    const result = [];

    for (let y = 0; y < u + v; ++y) {
        for (let x = -v; x < u; ++x) {
            if (indexedDB.cmp(reflect(x, y), [x, y])) {
                continue;
            }
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
