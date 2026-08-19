"""GenerateMultipointMacro_FromXML.py
Author, Hasreet Gill, Bassler Lab, HHMI/Princeton University
Code review/commenting with assistance from Claude (Sonnet 5, Anthropic)

Code generator that turns the multipoint-template NIS-Elements macro
XYCorrection_MultipointTimelapse_20X.mac (companion file, saved as
"xy20.mac" in the microscope's macros folder) into a ready-to-run macro
for however many wells/positions were actually defined that session.

The template macro is deliberately incomplete: it uses NUM_POINTS,
per-well bead-tracking arrays, per-well stage coordinates, and four
drift-correction code blocks that are never declared/filled in inside
the .mac file itself -- they're only placeholders, marked by exact
comment strings the macro author put there on purpose:
    "// INSERT_POINTS_HERE"                  (a line to be REPLACED)
    "// Object detection after shift"         (inserted AFTER this line)
    "// Compute median shift"                 (inserted AFTER this line)
    "// Update pre-shift"                     (inserted AFTER this line)
    "// First timepoint: just pre-shift"      (inserted AFTER this line)
plus the requirement that a "global double medianY" line and at least
one "#define" line exist to anchor the new global-array and NUM_POINTS
declarations. If any of those five marker strings, or the medianY/
#define anchors, are edited/removed in the .mac template, this script's
insertion points silently stop matching and that piece of generated code
just won't appear in the output -- so the template and this script must
be kept in sync.

This script reads "multipoints.xml" (exported from NIS-Elements' own
multipoint definition dialog -- see the checklist comment at the top of
the .mac template: points must be defined in Base, Surface, Base,
Surface, ... order, one pair per well), counts how many well-pairs that
implies, and generates exactly the pieces the template is missing:
NUM_POINTS, one pair of bead-position arrays per well, the literal
per-well stage coordinates, and one conditional code block per well for
each of the four drift-correction insertion points. It then writes out a
NEW, complete macro file, leaving the original template untouched so the
same template can be reused every time the well count changes.

Usage: run from the folder containing multipoints.xml and xy20.mac; it
writes xy20_points.mac alongside them, which is the file to actually
load and run in NIS-Elements.
"""

import xml.etree.ElementTree as ET

INDENT = "    " * 4  # Four tabs = 16 spaces -- matches the macro's nesting
                      # depth at the insertion points (inside the time loop,
                      # the multipoint loop, and the n>0/else block), so the
                      # generated if(m==i){...} blocks line up visually with
                      # the surrounding template code.

def parse_xml_positions(xml_file):
    """Reads the NIS-Elements multipoint XML export and splits its point
    list into Base (xB/yB/zB) and Surface (xS/yS/zS) coordinate lists.

    Points must have been defined in the microscope software in
    alternating Base, Surface, Base, Surface, ... order (one pair per
    well) -- this simply assigns even-indexed points to Base and
    odd-indexed points to Surface. The number of wells is len(xB)
    (equivalently len(xS)): half of however many <Point*> elements are
    present under the XML's "no_name" node.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    xB, yB, zB, xS, yS, zS = [], [], [], [], [], []
    points = root.find(".//no_name")
    point_list = [pt for pt in points if pt.tag.startswith("Point")]

    for i, pt in enumerate(point_list):
        x = float(pt.find("dXPosition").attrib["value"])
        y = float(pt.find("dYPosition").attrib["value"])
        z = float(pt.find("dZPosition").attrib["value"])
        if i % 2 == 0:
            xB.append(x)
            yB.append(y)
            zB.append(z)
        else:
            xS.append(x)
            yS.append(y)
            zS.append(z)
    return xB, yB, zB, xS, yS, zS

def indent_block(block):
    """Prefixes every non-blank line of BLOCK with INDENT, so a multi-line
    generated snippet (e.g. an if/{ /.../ } block) drops into the macro at
    the right nesting depth as plain text."""
    return "\n".join(f"{INDENT}{line}" if line else "" for line in block.splitlines())

def generate_macro_blocks(xB, yB, zB, xS, yS, zS):
    """Builds every piece of macro source code the template is missing,
    for exactly num = len(xB) wells:
      - define: the "#define NUM_POINTS <num>" line the template relies
        on but never declares itself.
      - globals: per-well pairs of bead-centroid arrays (CentreX/Y_1_MPi
        for the "before" snapshot, _2_MPi for the "after" snapshot) --
        one pair PER WELL is required (rather than one shared pair, as in
        the single-position acquisition macros) so each well's own
        previous bead snapshot survives, unclobbered by every other
        well's visits, across the outer time loop.
      - points: literal XpointsB/YpointsB/ZpointsB and XpointsS/YpointsS/
        ZpointsS[i] assignment statements -- the actual per-well stage
        coordinates read from the XML, replacing the macro's single
        "// INSERT_POINTS_HERE" placeholder line outright.
      - detect2 / shift / update1: one "if(m == i) { ... }" block per well
        for each of three drift-correction steps (post-shift bead
        detection, median-shift computation, and refreshing the "before"
        snapshot) -- the multipoint loop uses a single shared "m" index
        variable, so per-well behavior has to be dispatched with one
        if-block per well rather than array indexing directly by m
        (this macro dialect's arrays are declared with per-well NAMES,
        e.g. CentreX_1_MP3, not as a single array-of-arrays).
      - update1 is generated once and reused for BOTH the "Update
        pre-shift" marker (refreshing this well's snapshot after
        correcting this visit's drift) and the "First timepoint: just
        pre-shift" marker (taking the very first snapshot, before there's
        any previous one to correct against) -- both are literally the
        same "take a fresh bead snapshot for this well" operation.
    """
    num = len(xB)
    define = f"#define NUM_POINTS {num}\n"

    globals = "\n".join(
        f"global double CentreX_1_MP{i}[MAX_POINTS], CentreY_1_MP{i}[MAX_POINTS];\n"
        f"global double CentreX_2_MP{i}[MAX_POINTS], CentreY_2_MP{i}[MAX_POINTS];"
        for i in range(num)
    ) + "\n"

    points = ""
    for i in range(num):
        points += f"XpointsB[{i}] = {xB[i]:.4f};\nYpointsB[{i}] = {yB[i]:.4f};\nZpointsB[{i}] = {zB[i]:.4f};\n"
    for i in range(num):
        points += f"XpointsS[{i}] = {xS[i]:.4f};\nYpointsS[{i}] = {yS[i]:.4f};\nZpointsS[{i}] = {zS[i]:.4f};\n"

    detect2 = "\n\n".join(
        indent_block(f"""if(m == {i})
{{
    numPoints_2[{i}] = ObjectPositions(CentreX_2_MP{i}, CentreY_2_MP{i});
}}""")
        for i in range(num)
    )

    compute_shift = "\n\n".join(
        indent_block(f"""if(m == {i})
{{
    ComputeMedianShift(CentreX_1_MP{i}, CentreY_1_MP{i}, numPoints_1[{i}],
                       CentreX_2_MP{i}, CentreY_2_MP{i}, numPoints_2[{i}]);
}}""")
        for i in range(num)
    )

    update1 = "\n\n".join(
        indent_block(f"""if(m == {i})
{{
    numPoints_1[{i}] = ObjectPositions(CentreX_1_MP{i}, CentreY_1_MP{i});
}}""")
        for i in range(num)
    )

    return define, globals, points, detect2, compute_shift, update1

def insert_after_comment(lines, marker, content):
    """Finds the first line containing MARKER and inserts CONTENT as a new
    line immediately after it. Only the FIRST match is used, so each
    marker string must appear exactly once in the template at the spot
    meant to receive that particular generated block."""
    for i, line in enumerate(lines):
        if marker in line:
            lines.insert(i + 1, content + "\n")
            break

def modify_macro(xml_path, input_macro, output_macro):
    """Reads the multipoint XML + the template macro, generates every
    missing piece via generate_macro_blocks, splices it all into the
    template at the anchors described in this module's docstring, and
    writes the result to a NEW file (output_macro) -- the template
    (input_macro) itself is never modified, so it stays reusable."""
    xB, yB, zB, xS, yS, zS = parse_xml_positions(xml_path)
    define, globals, points, detect2, shift, update1 = generate_macro_blocks(xB, yB, zB, xS, yS, zS)

    with open(input_macro, 'r') as f:
        lines = f.readlines()

    # Insert #define NUM_POINTS after last #define
    idx_define = max(i for i, line in enumerate(lines) if "#define" in line)
    lines.insert(idx_define + 1, define + "\n")

    # Insert global arrays after known global section
    idx_global = next(i for i, line in enumerate(lines) if "global double medianY" in line)
    lines.insert(idx_global + 1, globals + "\n")

    # Replace // INSERT_POINTS_HERE
    idx_insert = next(i for i, line in enumerate(lines) if "// INSERT_POINTS_HERE" in line)
    lines[idx_insert] = points + "\n"

    # Insert code blocks below comments
    insert_after_comment(lines, "// Object detection after shift", detect2)
    insert_after_comment(lines, "// Compute median shift", shift)
    insert_after_comment(lines, "// Update pre-shift", update1)
    insert_after_comment(lines, "// First timepoint: just pre-shift", update1)

    with open(output_macro, 'w') as f:
        f.writelines(lines)

    print("✅ Macro updated and saved to:", output_macro)

# --- Run Script ---
modify_macro("multipoints.xml", "xy20.mac", "xy20_points.mac")
