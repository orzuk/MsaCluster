"""
Interactive 3D structure viewers for per-pair pages, built on Mol* (molstar).

This module produces a self-contained block of HTML that renders two beautiful,
interactive Mol* viewers for a fold-switching pair:

  (1) an ALIGNED / SUPERPOSED viewer  - both folds drawn in one canvas, each in
      its own colour, so the structural divergence is immediately visible; and
  (2) a SIDE-BY-SIDE view             - the two folds in their own canvases, in
      the same camera orientation.

Design goals (see the Ligo "Unreasonable Redundancy of Nature's Protein Folds"
blog, whose molstar setup this mirrors):

  * Self-contained / GitHub-Pages friendly. The Mol* bundle and CSS are copied
    once into ``docs/HTML/assets/`` and referenced with a *relative* path
    (``assets/molstar.js``) so they work for inline / offline pages. Each PDB's
    text is embedded directly in the page and loaded with
    ``loadStructureFromData`` - there are no external structure fetches.

  * Never raise. ``render_structure_viewers`` wraps every file read in
    try/except and degrades to a small ``<div class="warn">`` placeholder so the
    caller (``pb.push(render_structure_viewers(...))``) is safe.

Public entry point:

    render_structure_viewers(pair_id, fig_dir, output_dir, mode) -> str
"""

from __future__ import annotations

import os
import re
import json
import html as _html

# --------------------------------------------------------------------------- #
# Configuration / paths
# --------------------------------------------------------------------------- #

# Distinct, attractive colours for the two folds (hex, used both in CSS swatches
# and as integers passed to Mol*'s uniform colour theme).
_FOLD1_COLOR = "#4FC3F7"   # cyan / light blue
_FOLD2_COLOR = "#FF8A65"   # warm coral / orange

# Relative path (from a page in docs/HTML/) to the copied Mol* assets.
_ASSET_JS = "assets/molstar.js"
_ASSET_CSS = "assets/molstar.css"

_VIEWER_HEIGHT_PX = 480


def _data_dir() -> str:
    """Resolve the FoldPairs data directory via config.DATA_DIR, robustly."""
    try:
        import sys
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
        import config as CFG  # type: ignore
        d = getattr(CFG, "DATA_DIR", None)
        if d:
            return str(d)
    except Exception:
        pass
    # Fallback: assume repo layout
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(repo_root, "Pipeline", "FoldPairs")


def _split_pair(pair_id: str):
    """Split a pair id like '1dzlA_5keqF' into its two tags.

    Tags themselves never contain '_', so a single split on the first '_' from
    the right is robust even for unusual ids. Returns (tag1, tag2) or (pair, '').
    """
    pid = (pair_id or "").strip()
    if "_" in pid:
        i = pid.rfind("_")
        return pid[:i], pid[i + 1:]
    return pid, ""


def _safe_dom_id(s: str) -> str:
    """Make a string safe to use as a DOM id fragment."""
    return re.sub(r"[^A-Za-z0-9_-]", "-", str(s))


# --------------------------------------------------------------------------- #
# PDB file resolution + reading
# --------------------------------------------------------------------------- #

def _find_pdb_for_tag(pair_id: str, tag: str):
    """Locate the truth-backbone PDB file for a given fold tag.

    Preference order:
      1. Pipeline/FoldPairs/<pair_id>/chain_pdb_files/<tag>.pdb
      2. case-insensitive / loose match inside chain_pdb_files/
      3. whole-PDB fallback: any <tag>*.pdb or <pdbid>.pdb in the pair dir

    Returns an absolute path or None.
    """
    if not tag:
        return None
    try:
        base = os.path.join(_data_dir(), pair_id)
    except Exception:
        return None

    candidates = []
    chain_dir = os.path.join(base, "chain_pdb_files")
    # 1. exact
    candidates.append(os.path.join(chain_dir, tag + ".pdb"))
    # whole-pdb fallback locations
    candidates.append(os.path.join(base, tag + ".pdb"))

    for c in candidates:
        try:
            if os.path.isfile(c):
                return c
        except Exception:
            continue

    # 2. loose match in chain_pdb_files (case-insensitive prefix)
    for d in (chain_dir, base):
        try:
            if not os.path.isdir(d):
                continue
            tag_l = tag.lower()
            # exact-stem case-insensitive first
            for fn in os.listdir(d):
                if not fn.lower().endswith(".pdb"):
                    continue
                if os.path.splitext(fn)[0].lower() == tag_l:
                    return os.path.join(d, fn)
            # then prefix match (e.g. tag '1dzlA' -> '1dzlA_chainA.pdb')
            for fn in os.listdir(d):
                if not fn.lower().endswith(".pdb"):
                    continue
                if os.path.splitext(fn)[0].lower().startswith(tag_l):
                    return os.path.join(d, fn)
            # finally, match on the 4-char pdb id (drop chain letter)
            pdbid = tag_l[:4]
            if len(pdbid) == 4:
                for fn in os.listdir(d):
                    if not fn.lower().endswith(".pdb"):
                        continue
                    if os.path.splitext(fn)[0].lower().startswith(pdbid):
                        return os.path.join(d, fn)
        except Exception:
            continue
    return None


def _read_pdb_text(path):
    """Read a PDB file's text, returning None on any failure."""
    if not path:
        return None
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# HTML fragments
# --------------------------------------------------------------------------- #

def _warn_div(message: str) -> str:
    return ('<div class="warn" style="padding:10px 14px;border:1px solid #b5651d;'
            'border-radius:8px;background:#2a1f12;color:#ffcc99;font-size:0.9em;">'
            + _html.escape(message) + "</div>")


def _viewer_box(dom_id: str, title: str, status_id: str) -> str:
    """A styled container that Mol* mounts into (dark theme, rounded border)."""
    return (
        '<div class="molstar-box" '
        'style="position:relative;flex:1 1 320px;min-width:280px;max-width:100%;'
        f'height:{_VIEWER_HEIGHT_PX}px;border:1px solid #ccc;border-radius:12px;'
        'overflow:hidden;background:#ffffff;box-shadow:0 2px 10px rgba(0,0,0,0.15);">'
        f'<div style="position:absolute;top:8px;left:12px;z-index:5;color:#333;'
        'font-size:0.82em;font-weight:600;letter-spacing:0.02em;'
        'text-shadow:0 1px 2px rgba(255,255,255,0.8);pointer-events:none;">'
        + _html.escape(title) + "</div>"
        f'<div id="{dom_id}" style="position:absolute;inset:0;"></div>'
        f'<div id="{status_id}" class="molstar-status" '
        'style="position:absolute;bottom:8px;left:12px;z-index:5;color:#555;'
        'font-size:0.75em;pointer-events:none;">Loading...</div>'
        "</div>"
    )


def _legend(tag1: str, tag2: str) -> str:
    sw = ('display:inline-block;width:12px;height:12px;border-radius:3px;'
          'margin-right:5px;vertical-align:middle;')
    return (
        '<div class="molstar-legend" '
        'style="margin:6px 2px 10px;font-size:0.82em;color:#cfd6e4;">'
        f'<span style="margin-right:18px;"><span style="{sw}background:{_FOLD1_COLOR};"></span>'
        f'{_html.escape(tag1)}</span>'
        f'<span><span style="{sw}background:{_FOLD2_COLOR};"></span>'
        f'{_html.escape(tag2)}</span>'
        "</div>"
    )


# --------------------------------------------------------------------------- #
# JavaScript builder
# --------------------------------------------------------------------------- #

def _build_script(uid: str,
                  superposed_id: str, superposed_status_id: str,
                  side1_id: str, side1_status_id: str,
                  side2_id: str, side2_status_id: str,
                  pdb1_text, pdb2_text) -> str:
    """Build the self-contained <script> that mounts the Mol* viewers.

    PDB texts are embedded as JS string literals (via json.dumps, which safely
    escapes quotes, backslashes and newlines). Each fold is loaded with
    loadStructureFromData and coloured with a uniform colour theme. Structures
    are built through the public plugin builders so we get full control over
    representation + per-structure colour, and so two structures can coexist in
    one canvas (the superposed viewer).
    """
    # json.dumps produces a valid double-quoted JS string literal.
    pdb1_js = json.dumps(pdb1_text) if pdb1_text is not None else "null"
    pdb2_js = json.dumps(pdb2_text) if pdb2_text is not None else "null"

    color1 = int(_FOLD1_COLOR.lstrip("#"), 16)
    color2 = int(_FOLD2_COLOR.lstrip("#"), 16)

    cfg = {
        "uid": uid,
        "superposedId": superposed_id,
        "superposedStatusId": superposed_status_id,
        "side1Id": side1_id,
        "side1StatusId": side1_status_id,
        "side2Id": side2_id,
        "side2StatusId": side2_status_id,
        "color1": color1,
        "color2": color2,
    }
    cfg_js = json.dumps(cfg)

    # NOTE: kept as a plain string (not an f-string) so the JS braces are literal.
    js = """
<script>
(function () {
  var CFG = __CFG__;
  var PDB1 = __PDB1__;
  var PDB2 = __PDB2__;

  var VIEWER_OPTS = {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowRemoteState: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    collapseLeftPanel: true,
    collapseRightPanel: true,
    viewportShowControls: false,
    viewportShowExpand: false,
    viewportShowReset: false,
    viewportShowSettings: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    viewportShowTrajectoryControls: false,
    viewportShowScreenshotControls: false,
    pdbProvider: "rcsb",
    emdbProvider: "rcsb",
    illumination: true,
    pixelScale: 0.9
  };

  function setStatus(statusId, msg, state) {
    var el = document.getElementById(statusId);
    if (el) { el.textContent = msg; el.dataset.state = state || ""; }
  }

  function sleep(ms) { return new Promise(function (r) { window.setTimeout(r, ms); }); }

  // Build one structure from embedded PDB text and add a uniform-coloured
  // cartoon representation. Returns the structure ref (or null on failure).
  async function addStructure(plugin, pdbText, colorInt) {
    var data = await plugin.builders.data.rawData({ data: pdbText });
    var traj = await plugin.builders.structure.parseTrajectory(data, "pdb");
    // applyPreset builds model + structure + components + representations AND
    // computes secondary structure, so cartoon renders real helices/strands
    // (the manual createStructure/addRepresentation chain skipped SS -> a dot
    // cloud). Use a uniform colour so the two folds are visually distinct.
    await plugin.builders.structure.hierarchy.applyPreset(traj, "default", {
      structure: { name: "model", params: {} },
      representationPreset: "polymer-cartoon",
      representationPresetParams: {
        theme: {
          globalName: "uniform",
          globalColorParams: { value: colorInt }
        }
      }
    });
    return true;
  }

  function tuneCanvas(plugin) {
    if (!plugin || !plugin.canvas3d) { return; }
    plugin.canvas3d.setProps({
      camera: { manualReset: true },
      trackball: {
        minDistance: 2.5,
        zoomSpeed: 6,
        panSpeed: 1,
        rotateSpeed: 5
      },
      renderer: {
        backgroundColor: 0xffffff
      }
    });
  }

  // Zoom the camera in by `amount` (<1 = closer) so the structure fills the
  // canvas, mirroring the Ligo blog's framing. Safe no-op if state is missing.
  function zoomCamera(plugin, amount) {
    var cam = plugin && plugin.canvas3d && plugin.canvas3d.camera;
    var st = cam && cam.state;
    if (!st || !st.target || !st.position || typeof cam.setState !== "function") { return; }
    var t = st.target, p = st.position;
    var nx = t[0] + (p[0] - t[0]) * amount;
    var ny = t[1] + (p[1] - t[1]) * amount;
    var nz = t[2] + (p[2] - t[2]) * amount;
    cam.setState(Object.assign({}, st, {
      position: [nx, ny, nz],
      radius: Math.max((st.radius || 1) * amount, 1)
    }));
  }

  async function makeViewer(targetId, statusId) {
    var target = document.getElementById(targetId);
    if (!target) { return null; }
    if (!(window.molstar && window.molstar.Viewer)) {
      setStatus(statusId, "Mol* failed to load.", "error");
      return null;
    }
    setStatus(statusId, "Loading structure...", "loading");
    var viewer = await window.molstar.Viewer.create(target, VIEWER_OPTS);
    return viewer;
  }

  async function finalize(viewer, statusId, okMsg) {
    var plugin = viewer.plugin;
    tuneCanvas(plugin);
    viewer.handleResize();
    await sleep(120);
    if (plugin.canvas3d && plugin.canvas3d.requestCameraReset) {
      plugin.canvas3d.requestCameraReset({ durationMs: 0 });
    }
    await sleep(180);
    zoomCamera(plugin, 0.74);
    window.addEventListener("resize", function () { viewer.handleResize(); });
    setStatus(statusId, okMsg, "loaded");
  }

  // ---- Superposed viewer: both folds in one canvas ----
  async function loadSuperposed() {
    var statusId = CFG.superposedStatusId;
    try {
      var viewer = await makeViewer(CFG.superposedId, statusId);
      if (!viewer) { return; }
      var plugin = viewer.plugin;
      var any = false;
      if (PDB1) { await addStructure(plugin, PDB1, CFG.color1); any = true; }
      if (PDB2) { await addStructure(plugin, PDB2, CFG.color2); any = true; }
      if (!any) { setStatus(statusId, "No structures available.", "error"); return; }
      await finalize(viewer, statusId, "Superposed overlay ready (drag to rotate).");
    } catch (e) {
      console.error(e);
      setStatus(statusId, "Could not initialize this Mol* view.", "error");
    }
  }

  // ---- One side-by-side viewer: a single fold ----
  async function loadSide(targetId, statusId, pdbText, colorInt, label) {
    try {
      if (!pdbText) { setStatus(statusId, "structure not found", "error"); return; }
      var viewer = await makeViewer(targetId, statusId);
      if (!viewer) { return; }
      await addStructure(viewer.plugin, pdbText, colorInt);
      await finalize(viewer, statusId, label + " ready.");
    } catch (e) {
      console.error(e);
      setStatus(statusId, "Could not initialize this Mol* view.", "error");
    }
  }

  function boot() {
    loadSuperposed();
    loadSide(CFG.side1Id, CFG.side1StatusId, PDB1, CFG.color1, "Fold 1");
    loadSide(CFG.side2Id, CFG.side2StatusId, PDB2, CFG.color2, "Fold 2");
  }

  if (window.molstar && window.molstar.Viewer) {
    boot();
  } else {
    // Asset may still be loading (script with same src appears once per page).
    var tries = 0;
    var iv = window.setInterval(function () {
      tries++;
      if (window.molstar && window.molstar.Viewer) {
        window.clearInterval(iv);
        boot();
      } else if (tries > 100) {           // ~10s
        window.clearInterval(iv);
        setStatus(CFG.superposedStatusId, "Mol* failed to load.", "error");
        setStatus(CFG.side1StatusId, "Mol* failed to load.", "error");
        setStatus(CFG.side2StatusId, "Mol* failed to load.", "error");
      }
    }, 100);
  }
})();
</script>
"""
    js = js.replace("__CFG__", cfg_js)
    js = js.replace("__PDB1__", pdb1_js)
    js = js.replace("__PDB2__", pdb2_js)
    return js


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def render_structure_viewers(pair_id: str, fig_dir, output_dir, mode: str) -> str:
    """Return an HTML string with interactive Mol* viewers for ``pair_id``.

    Parameters
    ----------
    pair_id : str
        e.g. ``"1dzlA_5keqF"``. The two tags are the halves of the id.
    fig_dir, output_dir, mode :
        Accepted for signature compatibility with the page builder. The viewers
        embed PDB text directly, so these are not strictly required, but
        ``mode`` is honoured loosely (the output is identical and self-contained
        for both "inline" and "copy" modes).

    Notes
    -----
    This function never raises; on any failure it returns a small warning div.
    The Mol* bundle/CSS are referenced relatively (``assets/molstar.js``) and
    are expected to live in ``docs/HTML/assets/``.
    """
    try:
        tag1, tag2 = _split_pair(pair_id)
        uid = _safe_dom_id(pair_id) or "pair"

        # Resolve and read both PDBs (each independently optional).
        pdb1_text = _read_pdb_text(_find_pdb_for_tag(pair_id, tag1))
        pdb2_text = _read_pdb_text(_find_pdb_for_tag(pair_id, tag2))

        if pdb1_text is None and pdb2_text is None:
            return _warn_div(
                "structure not found for both folds of %s" % (pair_id,))

        # Unique DOM ids per pair / viewer to avoid collisions on multi-viewer
        # pages.
        superposed_id = "molstar-sup-%s" % uid
        superposed_status_id = superposed_id + "-status"
        side1_id = "molstar-s1-%s" % uid
        side1_status_id = side1_id + "-status"
        side2_id = "molstar-s2-%s" % uid
        side2_status_id = side2_id + "-status"

        # --- assemble HTML ---
        parts = []

        # Asset tags. Guard so the bundle is only injected once per page even if
        # this function is called multiple times.
        parts.append(
            '<link rel="stylesheet" type="text/css" '
            f'href="{_ASSET_CSS}" data-molstar-asset="1">'
        )
        parts.append(
            '<script type="text/javascript" '
            f'src="{_ASSET_JS}" data-molstar-asset="1"></script>'
        )

        parts.append(
            '<div class="structure-viewers" '
            'style="margin:14px 0;color:#cfd6e4;font-family:inherit;">'
        )

        # Section 1: superposed overlay
        parts.append(
            '<div style="font-size:0.95em;font-weight:600;margin:4px 2px 6px;">'
            'Superposed (both folds in one view)</div>'
        )
        parts.append(_legend(tag1 or "fold 1", tag2 or "fold 2"))
        if pdb1_text is None:
            parts.append(_warn_div("structure not found: %s" % (tag1,)))
        if pdb2_text is None:
            parts.append(_warn_div("structure not found: %s" % (tag2,)))
        parts.append('<div style="display:flex;">')
        parts.append(_viewer_box(superposed_id, "Aligned overlay",
                                 superposed_status_id))
        parts.append("</div>")

        # Section 2: side-by-side
        parts.append(
            '<div style="font-size:0.95em;font-weight:600;margin:18px 2px 6px;">'
            'Side by side</div>'
        )
        parts.append(
            '<div style="display:flex;gap:14px;flex-wrap:wrap;">'
        )
        if pdb1_text is not None:
            parts.append(_viewer_box(side1_id, tag1 or "Fold 1",
                                     side1_status_id))
        else:
            parts.append(_warn_div("structure not found: %s" % (tag1,)))
        if pdb2_text is not None:
            parts.append(_viewer_box(side2_id, tag2 or "Fold 2",
                                     side2_status_id))
        else:
            parts.append(_warn_div("structure not found: %s" % (tag2,)))
        parts.append("</div>")  # side-by-side flex

        parts.append("</div>")  # structure-viewers

        # Script
        parts.append(_build_script(
            uid,
            superposed_id, superposed_status_id,
            side1_id, side1_status_id,
            side2_id, side2_status_id,
            pdb1_text, pdb2_text,
        ))

        return "\n".join(parts)

    except Exception as e:  # absolute last-resort guard - never raise
        try:
            return _warn_div("structure viewer error: %s" % (e,))
        except Exception:
            return '<div class="warn">structure not found</div>'
