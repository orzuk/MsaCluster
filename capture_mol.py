#!/usr/bin/env python3
"""Capture one Mol* viewer element from a per-pair HTML page as a PNG.
Run under a virtual display:  xvfb-run -a python capture_mol.py <page.html> <element_id> <out.png> [W] [H]
Requires: pip install playwright ; playwright install chromium
          conda install -c conda-forge xorg-x11-server-xvfb mesalib   (for xvfb-run)
"""
import sys, os
from playwright.sync_api import sync_playwright

page_html = os.path.abspath(sys.argv[1])
element_id = sys.argv[2]
out_png = os.path.abspath(sys.argv[3])
W = int(sys.argv[4]) if len(sys.argv) > 4 else 1100
H = int(sys.argv[5]) if len(sys.argv) > 5 else 1000

GL_ARGS = [
    "--no-sandbox", "--use-gl=angle", "--use-angle=swiftshader",
    "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist",
    "--enable-webgl", "--disable-dev-shm-usage",
]

with sync_playwright() as p:
    # headless=False + a virtual X display (xvfb-run) is the reliable WebGL combo;
    # the swiftshader GL_ARGS act as a software-GL fallback.
    b = p.chromium.launch(headless=False, args=GL_ARGS)
    pg = b.new_page(viewport={"width": 1600, "height": 1200}, device_scale_factor=2)
    errs = []
    pg.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
    pg.goto("file://" + page_html, wait_until="load", timeout=120000)
    sel = "#" + element_id
    # enlarge + scroll the viewer box into view to trigger a full-res render
    pg.evaluate(
        """([id,w,h])=>{const e=document.getElementById(id);if(!e)return;
        const box=e.closest('.molstar-box')||e;box.style.width=w+'px';box.style.height=h+'px';
        box.scrollIntoView();window.dispatchEvent(new Event('resize'));}""",
        [element_id, W, H],
    )
    try:
        pg.wait_for_selector(sel + " canvas", timeout=90000)
    except Exception as e:
        print("NO CANVAS:", e)
        print("console errors:", " || ".join(errs[-6:]))
        b.close(); sys.exit(2)
    pg.wait_for_timeout(12000)                      # let Mol* load structures + superpose + settle
    pg.evaluate("()=>window.dispatchEvent(new Event('resize'))")
    pg.wait_for_timeout(3000)
    target = pg.eval_on_selector(sel, "e=>0")       # ensure present
    box = pg.query_selector(sel)
    # screenshot the bounding box (the .molstar-box) so we get the framed viewer
    handle = pg.evaluate_handle(
        "(id)=>document.getElementById(id).closest('.molstar-box')||document.getElementById(id)",
        element_id,
    )
    handle.as_element().screenshot(path=out_png)
    b.close()
    print("WROTE", out_png)
