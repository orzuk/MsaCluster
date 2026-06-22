#!/usr/bin/env python3
"""Capture a per-pair Mol* viewer (interactive 3D) as a static PNG for the paper.

The figure pipeline shows 3D via interactive Mol* viewers (fast, rendered in the
viewer's browser on their GPU). A PDF can't embed a live viewer, so for the few
pairs featured in the paper we screenshot the Mol* "aligned overlay" viewer
headlessly (Playwright + Chromium) under a virtual X display (xvfbwrapper) — no
GPU and no manual `xvfb-run` wrapper needed. Output matches the website exactly.

CLI:
    python utils/capture_molstar.py <pair_id>
    python utils/capture_molstar.py <page.html> <element_id> <out.png>

Requires: pip install playwright xvfbwrapper ; playwright install chromium
(The Xvfb binary must exist on the host; it does on the cluster. If a DISPLAY is
already set, or you wrap the call in `xvfb-run`, xvfbwrapper is not used.)
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_GL_ARGS = [
    "--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader",
    "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist",
    "--enable-webgl", "--disable-dev-shm-usage",
]


def default_paths(pair_id):
    """(per-pair HTML page, Mol* aligned-overlay element id, output PNG)."""
    html = os.path.join(REPO, "docs", "HTML", f"{pair_id}.html")
    element = f"molstar-sup-{pair_id}"
    out = os.path.join(REPO, "docs", "HTML", "figs", pair_id,
                       f"{pair_id}_two_structures_aligned.png")
    return html, element, out


def _clean_chrome(png_path):
    """White-out Mol* UI chrome (title, buttons, axis widget, status caption)."""
    from PIL import Image, ImageDraw
    im = Image.open(png_path).convert("RGB")
    w, h = im.size
    d = ImageDraw.Draw(im)
    white = (255, 255, 255)
    d.rectangle([0, 0, int(0.22 * w), int(0.055 * h)], fill=white)          # title (top-left)
    d.rectangle([int(0.86 * w), 0, w, int(0.10 * h)], fill=white)           # buttons (top-right)
    d.rectangle([0, int(0.90 * h), int(0.12 * w), h], fill=white)           # axis widget (bottom-left)
    d.rectangle([0, int(0.95 * h), int(0.60 * w), h], fill=white)           # "drag to rotate" caption
    im.save(png_path)


def capture_viewer(page_html, element_id, out_png, settle_ms=10000,
                   timeout_ms=180000, width=1100, height=1000):
    """Screenshot one Mol* viewer element from a per-pair HTML page to out_png."""
    from playwright.sync_api import sync_playwright
    page_html = os.path.abspath(page_html)
    out_png = os.path.abspath(out_png)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    sel = "#" + element_id
    with sync_playwright() as p:
        # headless=False under a virtual display is the reliable WebGL combo.
        b = p.chromium.launch(headless=False, args=_GL_ARGS)
        pg = b.new_page(viewport={"width": 1600, "height": 1200}, device_scale_factor=2)
        pg.goto("file://" + page_html, wait_until="load", timeout=120000)
        pg.evaluate(
            "([id,w,h])=>{const e=document.getElementById(id);if(!e)return;"
            "const box=e.closest('.molstar-box')||e;box.style.width=w+'px';box.style.height=h+'px';"
            "box.scrollIntoView();window.dispatchEvent(new Event('resize'));}",
            [element_id, width, height])
        pg.wait_for_selector(sel + " canvas", timeout=90000)
        # Hide every OTHER Mol* viewer so the software renderer only drives this one.
        pg.evaluate(
            "(id)=>{const keep=document.getElementById(id).closest('.molstar-box')||document.getElementById(id);"
            "document.querySelectorAll('.molstar-box').forEach(b=>{if(b!==keep)b.style.display='none';});}",
            element_id)
        pg.wait_for_timeout(settle_ms)           # let Mol* load structures + superpose + settle
        pg.evaluate("()=>window.dispatchEvent(new Event('resize'))")
        pg.wait_for_timeout(2000)
        pg.evaluate(
            "(id)=>{(document.getElementById(id).closest('.molstar-box')||document.getElementById(id))"
            ".scrollIntoView({block:'start'});}", element_id)
        pg.wait_for_timeout(1500)
        bb = pg.eval_on_selector(
            sel,
            "e=>{const b=(e.closest('.molstar-box')||e).getBoundingClientRect();"
            "return {x:Math.max(0,Math.round(b.left)),y:Math.max(0,Math.round(b.top)),"
            "w:Math.round(b.width),h:Math.round(b.height)};}")
        vw, vh = pg.viewport_size["width"], pg.viewport_size["height"]
        clip = {"x": bb["x"], "y": bb["y"],
                "width": min(bb["w"], vw - bb["x"]), "height": min(bb["h"], vh - bb["y"])}
        # A live WebGL canvas never goes "stable" -> element.screenshot() times out;
        # screenshot the page region (clip) instead.
        pg.screenshot(path=out_png, clip=clip, timeout=timeout_ms, animations="disabled")
        b.close()
    try:
        _clean_chrome(out_png)
    except Exception as e:
        print(f"[capture_molstar] chrome cleanup skipped: {e}")
    return out_png


def capture_pair(pair_id, **kw):
    """Capture the aligned-overlay Mol* viewer for a pair into its fig dir."""
    html, element, out = default_paths(pair_id)
    if not os.path.exists(html):
        raise FileNotFoundError(
            f"per-pair HTML not found: {html} (build it first, e.g. "
            f"gen_pair_html.py --pairs {pair_id} --make_figs)")
    return capture_viewer(html, element, out, **kw)


def run_under_display(fn, *a, **kw):
    """Call fn; if no X DISPLAY is set, start a virtual one (xvfbwrapper) first."""
    if os.environ.get("DISPLAY"):
        return fn(*a, **kw)
    try:
        from xvfbwrapper import Xvfb
    except ImportError:
        raise RuntimeError(
            "no DISPLAY and xvfbwrapper not installed; `pip install xvfbwrapper` "
            "or run under `xvfb-run`.")
    with Xvfb(width=1600, height=1200):
        return fn(*a, **kw)


def main(argv):
    if len(argv) == 1:                       # pair_id
        print("WROTE", run_under_display(capture_pair, argv[0]))
    elif len(argv) >= 3:                     # page.html element_id out.png
        print("WROTE", run_under_display(capture_viewer, argv[0], argv[1], argv[2]))
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main(sys.argv[1:])
