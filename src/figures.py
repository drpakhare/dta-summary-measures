"""Draws the two inline SVG figures for the deck from data/crp-dta.json and writes them
between <!--FIG:name--> ... <!--/FIG:name--> markers in slides/dta-summary-measures.html.
Run after src/make_data.py.
"""
import json, math, pathlib, re
from statistics import NormalDist

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = json.loads((ROOT / "data/crp-dta.json").read_text())
M = D["model"]; Z = NormalDist()

def pdf(x, mu, sd): return math.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))

def density_svg():
    W, H = 700, 400; L, R, T, B = 30, 20, 30, 70
    lo, hi = math.log(0.5), math.log(1000)
    X = lambda v: L + (v - lo) / (hi - lo) * (W - L - R)
    ymax = max(pdf(M["mu0"], M["mu0"], M["sd0"]), pdf(M["mu1"], M["mu1"], M["sd1"])) * 1.08
    # The no-SBI curve is scaled to the same peak area as SBI; both curves are densities.
    Y = lambda v: H - B - v / ymax * (H - T - B)
    xs = [lo + i * (hi - lo) / 240 for i in range(241)]
    def path(mu, sd):
        pts = " ".join(f"{X(x):.1f},{Y(pdf(x, mu, sd)):.1f}" for x in xs)
        return f"M{X(lo):.1f},{Y(0):.1f} L{pts} L{X(hi):.1f},{Y(0):.1f} Z"
    out = [f'<svg class="fig" viewBox="0 0 {W} {H}" width="{W}" role="img" aria-label="CRP distributions in infants with and without serious bacterial infection, with cut-offs at 20, 40 and 80 mg/L">']
    out.append(f'<path class="nd" d="{path(M["mu0"], M["sd0"])}"/>')
    out.append(f'<path class="dz" d="{path(M["mu1"], M["sd1"])}"/>')
    out.append(f'<line class="axis" x1="{L}" x2="{W-R}" y1="{Y(0)}" y2="{Y(0)}"/>')
    for v in (1, 5, 10, 20, 40, 80, 200, 500):
        out.append(f'<text class="tick" x="{X(math.log(v)):.1f}" y="{Y(0)+26:.1f}" text-anchor="middle">{v}</text>')
    out.append(f'<text class="tick" x="{(W)/2}" y="{H-12}" text-anchor="middle">CRP, mg/L (log scale)</text>')
    for t in (20, 40, 80):
        x = X(math.log(t))
        out.append(f'<line class="cut" x1="{x:.1f}" x2="{x:.1f}" y1="{T+10}" y2="{Y(0)}"/>')
        out.append(f'<text class="cutl" x="{x:.1f}" y="{T}" text-anchor="middle">{t}</text>')
    out.append(f'<text class="lbl" x="{X(math.log(2.2)):.1f}" y="{Y(pdf(M["mu0"],M["mu0"],M["sd0"]))-10:.1f}" text-anchor="middle">No SBI</text>')
    out.append(f'<text class="lbl dzl" x="{X(math.log(200)):.1f}" y="{Y(pdf(math.log(200),M["mu1"],M["sd1"]))-14:.1f}" text-anchor="middle">SBI</text>')
    out.append('</svg>')
    return "\n".join(out)

def roc_svg():
    W, H = 470, 470; L, R, T, B = 80, 20, 20, 70
    X = lambda f: L + f * (W - L - R)
    Y = lambda s: H - B - s * (H - T - B)
    out = [f'<svg class="fig" viewBox="0 0 {W} {H}" width="{W}" role="img" aria-label="The five studies plotted in ROC space">']
    for g in (0.25, 0.5, 0.75, 1):
        out.append(f'<line class="grid" x1="{X(0)}" x2="{X(1)}" y1="{Y(g):.1f}" y2="{Y(g):.1f}"/>')
        out.append(f'<line class="grid" x1="{X(g):.1f}" x2="{X(g):.1f}" y1="{Y(0)}" y2="{Y(1)}"/>')
    out.append(f'<line class="axis" x1="{X(0)}" x2="{X(1)}" y1="{Y(0)}" y2="{Y(0)}"/>')
    out.append(f'<line class="axis" x1="{X(0)}" x2="{X(0)}" y1="{Y(0)}" y2="{Y(1)}"/>')
    out.append(f'<line class="chance" x1="{X(0)}" x2="{X(1)}" y1="{Y(0)}" y2="{Y(1)}"/>')
    for v in (0, 0.5, 1):
        out.append(f'<text class="tick" x="{X(v):.1f}" y="{Y(0)+26}" text-anchor="middle">{v:g}</text>')
        out.append(f'<text class="tick" x="{X(0)-10}" y="{Y(v)+7:.1f}" text-anchor="end">{v:g}</text>')
    out.append(f'<text class="tick" x="{(X(0)+X(1))/2}" y="{H-14}" text-anchor="middle">1 − specificity</text>')
    out.append(f'<text class="tick" x="22" y="{(Y(0)+Y(1))/2}" text-anchor="middle" transform="rotate(-90 22 {(Y(0)+Y(1))/2})">Sensitivity</text>')
    pts = " ".join(f'{X(c["fpr"]):.1f},{Y(c["tpr"]):.1f}' for c in D["curve"]) if "curve" in D else ""
    # curve from model
    cur = []
    for i in range(-60, 200):
        t = math.exp(i / 25)
        f = 1 - Z.cdf((math.log(t) - M["mu0"]) / M["sd0"]); s = 1 - Z.cdf((math.log(t) - M["mu1"]) / M["sd1"])
        cur.append(f"{X(f):.1f},{Y(s):.1f}")
    out.append(f'<polyline class="curve" points="{" ".join(cur)}"/>')
    cls = {20: "t20", 40: "t40", 80: "t80"}
    for s in D["studies"]:
        sn = s["TP"] / (s["TP"] + s["FN"]); fpr = s["FP"] / (s["FP"] + s["TN"])
        r = 6 + 8 * math.sqrt(s["N"] / 1100)
        out.append(f'<circle class="pt {cls[s["threshold"]]}" cx="{X(fpr):.1f}" cy="{Y(sn):.1f}" r="{r:.1f}"/>')
        if s["id"] == "C":
            out.append(f'<text class="pl" x="{X(fpr):.1f}" y="{Y(sn)-r-6:.1f}" text-anchor="middle">{s["id"]}</text>')
        else:
            out.append(f'<text class="pl" x="{X(fpr)+r+5:.1f}" y="{Y(sn)+7:.1f}">{s["id"]}</text>')
    # naive separate pooling (fixed-effect logit), for illustration
    def pool(pairs):
        num = den = 0
        for x, n in pairs:
            v = 1 / x + 1 / (n - x); num += math.log(x / (n - x)) / v; den += 1 / v
        return 1 / (1 + math.exp(-num / den))
    psn = pool([(s["TP"], s["TP"] + s["FN"]) for s in D["studies"]])
    psp = pool([(s["TN"], s["TN"] + s["FP"]) for s in D["studies"]])
    x, y = X(1 - psp), Y(psn)
    out.append(f'<path class="naive" d="M{x-10:.1f},{y-10:.1f} L{x+10:.1f},{y+10:.1f} M{x-10:.1f},{y+10:.1f} L{x+10:.1f},{y-10:.1f}"/>')
    out.append(f'<text class="pl naivel" x="{x+16:.1f}" y="{y+26:.1f}">pooled separately</text>')
    out.append('</svg>')
    return "\n".join(out), psn, psp

deck = ROOT / "slides/dta-summary-measures.html"
h = deck.read_text()
roc, psn, psp = roc_svg()
for name, svg in (("density", density_svg()), ("roc", roc)):
    h = re.sub(rf"<!--FIG:{name}-->.*?<!--/FIG:{name}-->", lambda m: f"<!--FIG:{name}-->\n{svg}\n<!--/FIG:{name}-->", h, flags=re.S)
deck.write_text(h)
print(f"naive pooled Sn {psn:.3f} Sp {psp:.3f}")
