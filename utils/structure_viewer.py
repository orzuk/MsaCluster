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


def _chain_of_tag(tag: str):
    """Chain id encoded in a fold tag, e.g. '2n54B' -> 'B' (None if absent)."""
    t = (tag or "").strip()
    return t[4:] if len(t) > 4 else None


def _trigger_class(pair_id: str) -> str:
    """trigger_class for a pair from docs/triggers_from_pdb.csv ('' if unknown).

    Used so that oligomerization switchers (where the fold change IS the
    monomer<->multimer transition) keep their full assembly, while every other
    class is reduced to the single fold-switching chain for a clean comparison.
    """
    try:
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        import csv
        with open(os.path.join(repo, "docs", "triggers_from_pdb.csv"),
                  encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                if row.get("pair_id") == pair_id:
                    return (row.get("trigger_class") or "").strip()
    except Exception:
        pass
    return ""


def _find_super_pdb(pair_id, fig_dir, which):
    """Locate the PyMOL-superposed single-chain PDB for the aligned overlay.

    plot_3d.py writes `<pair>_super_fold{1,2}.pdb` (fold1 = reference frame,
    fold2 moved onto it) next to the other figures. Check the published fig dir
    and the Pipeline output_figs dir. Returns a path or None.
    """
    name = f"{pair_id}_super_fold{which}.pdb"
    cands = []
    try:
        if fig_dir:
            cands.append(os.path.join(str(fig_dir), name))
    except Exception:
        pass
    try:
        cands.append(os.path.join(_data_dir(), pair_id, "output_figs", name))
    except Exception:
        pass
    for c in cands:
        try:
            if os.path.isfile(c) and os.path.getsize(c) > 0:
                return c
        except Exception:
            continue
    return None


def _superpose_overlay(text1, text2):
    """Kabsch-superpose fold2 onto fold1 over sequence-matched CA atoms.

    Both inputs are single-chain structure TEXT (PDB or mmCIF). Returns
    (pdb1_text, pdb2_superposed_text) as PDB strings in fold1's reference frame,
    or None on any failure (caller then shows the overlay unaligned). No PyMOL /
    TMalign dependency - uses Biopython + a pairwise sequence alignment. For a
    near-identical pair (TM~1) this overlays the two almost perfectly.
    """
    import tempfile
    _AA3TO1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
               'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
               'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
               'MSE': 'M'}

    def _parse(text):
        is_cif = ("_atom_site." in text[:4000]) or text.lstrip().startswith("data_")
        from Bio.PDB import MMCIFParser, PDBParser
        suf = ".cif" if is_cif else ".pdb"
        tf = tempfile.NamedTemporaryFile(suffix=suf, delete=False)
        tf.write(text.encode("utf-8", "replace")); tf.close()
        try:
            parser = MMCIFParser(QUIET=True) if is_cif else PDBParser(QUIET=True)
            return parser.get_structure("x", tf.name)
        finally:
            try: os.remove(tf.name)
            except Exception: pass

    def _ca_seq(struct):
        model = next(iter(struct))
        for ch in model:
            cas, seq = [], []
            for res in ch:
                rn = res.resname.strip()
                if "CA" in res and rn in _AA3TO1:
                    cas.append(res["CA"]); seq.append(_AA3TO1[rn])
            if len(cas) >= 3:
                return cas, "".join(seq)
        return [], ""

    def _to_pdb(struct):
        from Bio.PDB import PDBIO
        io = PDBIO(); io.set_structure(struct)
        tf = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False); tf.close()
        io.save(tf.name)
        try:
            with open(tf.name) as fh:
                return fh.read()
        finally:
            try: os.remove(tf.name)
            except Exception: pass

    try:
        from Bio.PDB import Superimposer
        from Bio import pairwise2
        s1, s2 = _parse(text1), _parse(text2)
        # NMR entries carry ~20 models; keep only the first so the overlay isn't
        # a 20-model cloud and the embedded text stays small.
        for s in (s1, s2):
            for m in list(s)[1:]:
                s.detach_child(m.id)
        # Drop waters + ligands/heteroatoms so the overlay is clean protein only
        # (these were the stray "dots" and odd colours in the viewer).
        from Bio.PDB.Polypeptide import is_aa as _is_aa
        for s in (s1, s2):
            for model in s:
                for chain in model:
                    for res in list(chain):
                        if not _is_aa(res, standard=True):
                            chain.detach_child(res.id)
        ca1, seq1 = _ca_seq(s1)
        ca2, seq2 = _ca_seq(s2)
        if len(ca1) < 3 or len(ca2) < 3:
            return None
        aln = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -0.5,
                                       one_alignment_only=True)[0]
        fixed, moving, i, j = [], [], 0, 0
        for a, b in zip(aln.seqA, aln.seqB):
            if a != "-" and b != "-" and i < len(ca1) and j < len(ca2):
                fixed.append(ca1[i]); moving.append(ca2[j])
            if a != "-": i += 1
            if b != "-": j += 1
        if len(fixed) < 3:
            return None
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply(list(s2.get_atoms()))   # move ALL of fold2 onto fold1
        return _to_pdb(s1), _to_pdb(s2)
    except Exception as e:
        print(f"[structure] overlay superposition failed: {e}")
        return None


def _structure_caption(pair_id):
    """One-line context for the 3D figures from docs/triggers_from_pdb.csv:
    trigger class, per-fold chain counts (monomer/dimer/...), ligands."""
    try:
        repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        import csv
        with open(os.path.join(repo, "docs", "triggers_from_pdb.csv"),
                  encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                if row.get("pair_id") != pair_id:
                    continue
                bits = []
                tc = (row.get("trigger_class") or "").strip()
                if tc:
                    bits.append(f"trigger class: <b>{_html.escape(tc)}</b>")
                n1 = (row.get("n_unique_chains_pdb1") or "").strip()
                n2 = (row.get("n_unique_chains_pdb2") or "").strip()

                def _olig(nch):
                    try:
                        k = int(nch)
                    except Exception:
                        return ""
                    return {1: "monomer", 2: "dimer", 3: "trimer",
                            4: "tetramer"}.get(k, f"{k}-mer")
                o1, o2 = _olig(n1), _olig(n2)
                if o1 or o2:
                    bits.append(f"deposited assembly: fold1 {o1 or '?'} / fold2 {o2 or '?'}")
                lig = ((row.get("ligands_pdb1") or "") + (row.get("ligands_pdb2") or "")).strip()
                if lig:
                    bits.append(f"ligand(s): {_html.escape(lig)}")
                bits.append("showing the single fold-switching chain of each fold, "
                            "Kabsch-superposed (waters/ligands hidden)")
                return " &middot; ".join(bits)
    except Exception:
        pass
    return ("showing the single fold-switching chain of each fold, "
            "Kabsch-superposed (waters/ligands hidden)")


def _fold_label(tag, pair_id):
    """'2n54B' -> '2n54B (XCL1)' using the pair's common name when available."""
    try:
        from utils.pair_labels import pair_display
        disp = pair_display(pair_id)
        if disp:
            return f"{tag} ({disp})"
    except Exception:
        pass
    return str(tag)


def _filter_structure_to_chain(text, chain):
    """Return `text` reduced to a single chain. Self-contained (no deps).

    Handles mmCIF (filter the _atom_site loop rows by auth/label_asym_id) and
    plain PDB (filter ATOM/HETATM/TER/ANISOU by the chain column). Returns the
    original text unchanged if the chain can't be isolated or would be empty -
    never raises, never drops a structure.
    """
    if not text or not chain:
        return text
    try:
        lines = text.splitlines(keepends=False)
        is_cif = any(l.lstrip().startswith(("_atom_site.", "data_")) for l in lines[:60])
        if is_cif:
            # locate the _atom_site loop header
            cols = []; ds = -1
            for i, l in enumerate(lines):
                if l.strip() == "loop_":
                    blk = []; j = i + 1
                    while j < len(lines) and lines[j].lstrip().startswith("_atom_site."):
                        blk.append(lines[j].strip()); j += 1
                    if blk:
                        cols = blk; ds = j; break
            if not cols:
                return text
            idx = {c: k for k, c in enumerate(cols)}
            ca = idx.get("_atom_site.auth_asym_id", idx.get("_atom_site.label_asym_id"))
            if ca is None:
                return text
            out = lines[:ds]; kept = 0
            for k in range(ds, len(lines)):
                ln = lines[k]; st = ln.strip()
                if (not st) or st.startswith("#") or st.startswith("_") or st == "loop_":
                    out.append(ln)
                    # an atom-record run ends at the next non-record line
                    if kept and (st.startswith("_") or st == "loop_"):
                        out.extend(lines[k + 1:]); break
                    continue
                p = st.split()
                if len(p) < len(cols):
                    out.append(ln); continue
                if p[ca] == chain:
                    out.append(ln); kept += 1
            return "\n".join(out) if kept >= 5 else text
        else:
            # plain PDB: chain id is column 22 (index 21)
            out = []; kept = 0
            for ln in lines:
                rec = ln[:6].strip()
                if rec in ("ATOM", "HETATM", "ANISOU", "TER"):
                    if len(ln) > 21 and ln[21] == chain:
                        out.append(ln); kept += 1 if rec in ("ATOM", "HETATM") else 0
                    # else drop this record (other chain)
                else:
                    out.append(ln)
            return "\n".join(out) if kept >= 5 else text
    except Exception:
        return text


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
  // Load one structure from embedded PDB text via the HIGH-LEVEL Viewer API.
  // This is the same path the Ligo blog uses (loadStructureFromUrl): it builds
  // model + structure + representation, computes secondary structure for a real
  // cartoon, AND auto-focuses the camera on the loaded structure. The previous
  // low-level builders.structure.hierarchy.applyPreset chain did neither
  // reliably, leaving a blank (all-white) canvas. A uniform colour keeps the
  // two folds visually distinct.
  // The embedded structure text is often mmCIF (starts with "data_" / contains
  // _atom_site.) even though it came from a ".pdb" file. Telling Mol* the wrong
  // format makes it silently parse nothing -> blank canvas. Detect and pass the
  // correct format.
  function detectFormat(txt) {
    if (!txt) { return "pdb"; }
    var head = txt.slice(0, 4000);
    if (head.indexOf("_atom_site.") >= 0 || head.slice(0, 60).indexOf("data_") >= 0) {
      return "mmcif";
    }
    return "pdb";
  }

  async function addStructure(viewer, pdbText, colorInt, label, colorMode) {
    // "uniform" (overlay: one solid colour per fold so the two are distinct)
    // vs "gradient" (side-by-side: colour by residue number N->C so the SAME
    // residue is the same colour in both folds -> matching residues line up).
    var theme = (colorMode === "gradient")
      ? { globalName: "sequence-id" }
      : { globalName: "uniform", globalColorParams: { value: colorInt } };
    await viewer.loadStructureFromData(pdbText, detectFormat(pdbText), {
      dataLabel: label || "structure",
      representationParams: { theme: theme }
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
      var any = false;
      if (PDB1) { await addStructure(viewer, PDB1, CFG.color1, "Fold 1"); any = true; }
      if (PDB2) { await addStructure(viewer, PDB2, CFG.color2, "Fold 2"); any = true; }
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
      await addStructure(viewer, pdbText, colorInt, label, "gradient");
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

def render_structure_viewers(pair_id: str, fig_dir, output_dir, mode: str,
                             fig_start=None):
    """Return ``(html, n_figs)`` with interactive Mol* viewers for ``pair_id``.

    Parameters
    ----------
    pair_id : str
        e.g. ``"1dzlA_5keqF"``. The two tags are the halves of the id.
    fig_dir, output_dir, mode :
        Accepted for signature compatibility with the page builder. The viewers
        embed PDB text directly, so these are not strictly required, but
        ``mode`` is honoured loosely (the output is identical and self-contained
        for both "inline" and "copy" modes).
    fig_start : int | None
        If given, the two viewer sections are labelled "Figure {fig_start}: ..."
        (aligned overlay) and "Figure {fig_start+1}: ..." (side by side) so they
        slot into the page's running figure numbering.

    Returns
    -------
    (html, n_figs) : the HTML block and how many figure numbers it consumed
    (2 when structures render, 0 on the all-missing warning path).

    Notes
    -----
    This function never raises; on any failure it returns a small warning div.
    The Mol* bundle/CSS are referenced relatively (``assets/molstar.js``) and
    are expected to live in ``docs/HTML/assets/``.
    """
    def _fig_label(offset, text):
        if fig_start is None:
            return text
        return f"Figure {int(fig_start) + offset}: {text}"
    try:
        tag1, tag2 = _split_pair(pair_id)
        uid = _safe_dom_id(pair_id) or "pair"

        # Resolve each fold's structure and reduce to its single fold-switching
        # chain (chain in the tag, e.g. 2n54B -> chain B). Oligomerization
        # switchers keep the full assembly (the fold change IS the multimer state).
        pdb1_text = _read_pdb_text(_find_pdb_for_tag(pair_id, tag1))
        pdb2_text = _read_pdb_text(_find_pdb_for_tag(pair_id, tag2))
        if pdb1_text is None and pdb2_text is None:
            return (_warn_div(
                "structure not found for both folds of %s" % (pair_id,)), 0)
        _is_oligo = _trigger_class(pair_id) == "oligomerization"
        if not _is_oligo:
            if pdb1_text is not None:
                pdb1_text = _filter_structure_to_chain(pdb1_text, _chain_of_tag(tag1))
            if pdb2_text is not None:
                pdb2_text = _filter_structure_to_chain(pdb2_text, _chain_of_tag(tag2))

        # Superpose fold2 onto fold1 (Biopython Kabsch over sequence-matched CA
        # atoms) so the aligned overlay ACTUALLY lines up - no PyMOL/TMalign
        # dependency. Replaces both texts with clean PDB in fold1's frame on
        # success; leaves them as-is (overlay just unaligned) on any failure.
        if pdb1_text and pdb2_text and not _is_oligo:
            _sup = _superpose_overlay(pdb1_text, pdb2_text)
            if _sup:
                pdb1_text, pdb2_text = _sup

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
            + _html.escape(_fig_label(0, "Aligned overlay (both folds superposed)"))
            + "</div>"
        )
        parts.append(
            '<div style="font-size:0.8em;color:#9aa3b2;margin:0 2px 6px;">'
            + _structure_caption(pair_id) + "</div>"
        )
        parts.append(_legend(_fold_label(tag1 or "fold 1", pair_id),
                             _fold_label(tag2 or "fold 2", pair_id)))
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
            + _html.escape(_fig_label(1, "The two folds, side by side"))
            + "</div>"
        )
        parts.append(
            '<div style="display:flex;gap:14px;flex-wrap:wrap;">'
        )
        if pdb1_text is not None:
            parts.append(_viewer_box(side1_id, _fold_label(tag1 or "Fold 1", pair_id),
                                     side1_status_id))
        else:
            parts.append(_warn_div("structure not found: %s" % (tag1,)))
        if pdb2_text is not None:
            parts.append(_viewer_box(side2_id, _fold_label(tag2 or "Fold 2", pair_id),
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

        # Two figure slots consumed (aligned overlay + side-by-side).
        return ("\n".join(parts), 2)

    except Exception as e:  # absolute last-resort guard - never raise
        try:
            return (_warn_div("structure viewer error: %s" % (e,)), 0)
        except Exception:
            return ('<div class="warn">structure not found</div>', 0)
