// Shared arithmetic for the labs. Wilson score intervals for proportions;
// log-method intervals for likelihood ratios and the diagnostic odds ratio.
// If any cell is zero, 0.5 is added to every cell for LRs and DOR only.
(function () {
  var Z95 = 1.959964;
  function wilson(x, n) {
    if (!n) return [NaN, NaN];
    var p = x / n, z = Z95, den = 1 + z * z / n;
    var c = (p + z * z / (2 * n)) / den;
    var h = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den;
    return [Math.max(0, c - h), Math.min(1, c + h)];
  }
  function measures(tp, fp, fn, tn) {
    var d = tp + fn, nd = fp + tn, pos = tp + fp, neg = fn + tn, n = d + nd;
    var r = { n: n, d: d, nd: nd, prev: d / n };
    r.sn = tp / d; r.snCI = wilson(tp, d);
    r.sp = tn / nd; r.spCI = wilson(tn, nd);
    r.ppv = tp / pos; r.ppvCI = wilson(tp, pos);
    r.npv = tn / neg; r.npvCI = wilson(tn, neg);
    var zero = (tp === 0 || fp === 0 || fn === 0 || tn === 0);
    var a = tp, b = fp, c = fn, e = tn;
    if (zero) { a += .5; b += .5; c += .5; e += .5; }
    var sn = a / (a + c), sp = e / (e + b);
    r.corrected = zero;
    r.lrp = sn / (1 - sp);
    var s1 = Math.sqrt(1 / a - 1 / (a + c) + 1 / b - 1 / (b + e));
    r.lrpCI = [r.lrp * Math.exp(-Z95 * s1), r.lrp * Math.exp(Z95 * s1)];
    r.lrn = (1 - sn) / sp;
    var s2 = Math.sqrt(1 / c - 1 / (a + c) + 1 / e - 1 / (b + e));
    r.lrnCI = [r.lrn * Math.exp(-Z95 * s2), r.lrn * Math.exp(Z95 * s2)];
    r.dor = (a * e) / (b * c);
    var s3 = Math.sqrt(1 / a + 1 / b + 1 / c + 1 / e);
    r.dorCI = [r.dor * Math.exp(-Z95 * s3), r.dor * Math.exp(Z95 * s3)];
    return r;
  }
  function post(pre, lr) { var o = pre / (1 - pre) * lr; return o / (1 + o); }
  function pct(x, dp) { return isFinite(x) ? (100 * x).toFixed(dp == null ? 1 : dp) + '%' : 'n/a'; }
  function num(x, dp) {
    if (!isFinite(x)) return 'n/a';
    if (dp == null) dp = x >= 100 ? 0 : x >= 10 ? 1 : 2;
    return x.toFixed(dp);
  }
  function ci(lo, hi, f) { return isFinite(lo) ? f(lo) + ' to ' + f(hi) : ''; }
  // Normal CDF (Abramowitz and Stegun 7.1.26 via erf), accurate to about 1e-7
  function phi(x) {
    var t = 1 / (1 + 0.3275911 * Math.abs(x) / Math.SQRT2);
    var y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x / 2);
    return x >= 0 ? (1 + y) / 2 : (1 - y) / 2;
  }
  window.DTA = { wilson: wilson, measures: measures, post: post, pct: pct, num: num, ci: ci, phi: phi };
})();
