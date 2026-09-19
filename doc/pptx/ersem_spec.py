"""
ERSEM schematic specification (logic layer), traced node-by-node and
edge-by-edge from docs/images/ERSEM.png (original ERSEM, Butenschoen 2016).

This file is the SINGLE SOURCE OF TRUTH for the diagram logic.  The pptx
builder (build_ersem_v2.py) only renders what is declared here, so the
connection logic can be reviewed and corrected independently of styling.

Coordinate system: centimetres on a 33.0 x 30.0 cm slide, origin top-left.
Each node has a centre (cx, cy) and a size (w, h).  Ports are derived
automatically (left/right/top/bottom mid-points) by the builder.

Edge colours (verified against the original figure):
    red    inorganic nutrient exchange (uptake / remineralisation)
    blue   carbonate system / DIC / O2 / gas exchange
    green  organic matter flow (DOM / POM production and transfer)
    black  predation / grazing (trophic transfer)
    grey   sinking flux (pelagic -> benthic)
"""

# ---------------------------------------------------------------------------
# NODE KINDS -> default styling handled by builder
#   capsule  : white/grey rounded inorganic species (Fe, Si, ...)
#   spec      : carbonate speciation dark box (DIC species)
#   gas       : small white oval (pCO2, DIC, O2, pH/Omega, TA)
#   phyto     : green ellipse (phytoplankton species)
#   bact      : orange ellipse (bacteria)
#   zoo       : blue-grey ellipse (zooplankton / zoobenthos)
#   hetero    : teal ellipse (heterotrophic flagellates)
#   pool      : light grey ellipse (DOM / Particulates)
#   group     : translucent rounded container (background)
# ---------------------------------------------------------------------------

# Band geometry --------------------------------------------------------------
BANDS = {
    "atmos":   dict(y0=0.6,  y1=2.0,  label="Atmosphere"),
    "pelagic": dict(y0=2.0,  y1=16.0, label="Pelagic"),
    "benthic": dict(y0=16.0, y1=29.4, label="Benthic"),
}
SEP_X = 9.9          # dashed Inorganics|Organics separator
HEADERS = {"Inorganics": 5.0, "Organics": 20.0}

# Group containers (drawn first, behind nodes) ------------------------------
GROUPS = [
    dict(name="Phytoplankton", x=11.3, y=2.9,  w=9.4, h=2.9, label_pos="tl"),
    dict(name="Microbes_pel",  x=11.3, y=6.4,  w=3.4, h=2.1, label="Microbes", label_pos="tl"),
    dict(name="Consumers",     x=11.3, y=9.1,  w=9.4, h=2.4, label="Consumers", label_pos="br"),
    dict(name="Zoobenthos",    x=11.3, y=17.0, w=6.6, h=3.8, label="Zoobenthos", label_pos="tr"),
    dict(name="Microbes_ben",  x=11.3, y=21.2, w=3.6, h=3.6, label="Microbes", label_pos="br"),
]

# ---------------------------------------------------------------------------
# NODES
# ---------------------------------------------------------------------------
NODES = {
    # ---- PELAGIC carbonate system ----
    "pco2_atm_p": dict(cx=2.2, cy=1.3, w=1.7, h=1.1, kind="gas", text="pCO₂", band="atmos"),
    "pco2_p":     dict(cx=2.2, cy=3.5, w=1.8, h=1.4, kind="gas", text="pCO₂"),
    "dic_p":      dict(cx=4.4, cy=3.2, w=1.5, h=1.2, kind="gas", text="DIC"),
    "o2_p":       dict(cx=4.4, cy=4.9, w=1.5, h=1.2, kind="gas", text="O₂"),
    "dicspec_p":  dict(cx=2.2, cy=6.6, w=2.2, h=2.6, kind="spec",
                       text="H₂CO₃\n⇅\nHCO₃⁻\n⇅\nCO₃²⁻"),
    "ph_p":       dict(cx=4.4, cy=6.6, w=1.6, h=1.3, kind="gas", text="pH\nΩ"),
    "ta_p":       dict(cx=2.2, cy=9.7, w=1.6, h=1.3, kind="gas", text="TA"),

    # ---- PELAGIC nutrients ----
    "Fe":  dict(cx=6.4, cy=4.6, w=1.6, h=1.15, kind="capsule", text="Fe"),
    "Si":  dict(cx=6.4, cy=5.9, w=1.6, h=1.15, kind="capsule", text="Si"),
    "PO4": dict(cx=6.4, cy=7.2, w=1.6, h=1.15, kind="capsule", text="PO₄"),
    "NO3": dict(cx=6.4, cy=8.5, w=1.6, h=1.15, kind="capsule", text="NO₃"),
    "NH4": dict(cx=6.4, cy=9.8, w=1.6, h=1.15, kind="capsule", text="NH₄"),

    # ---- PELAGIC phytoplankton (4 PFTs) ----
    "Picoph":  dict(cx=12.85, cy=4.35, w=2.0, h=1.5, kind="phyto", text="Picoph."),
    "Nanoph":  dict(cx=15.05, cy=4.35, w=2.0, h=1.5, kind="phyto", text="Nanoph."),
    "Microph": dict(cx=17.25, cy=4.35, w=2.0, h=1.5, kind="phyto", text="Microph."),
    "Diatoms": dict(cx=19.45, cy=4.35, w=2.0, h=1.5, kind="phyto", text="Diatoms"),

    # ---- PELAGIC microbes ----
    "Bacteria_p": dict(cx=13.0, cy=7.55, w=2.4, h=1.2, kind="bact", text="Bacteria"),

    # ---- PELAGIC consumers ----
    "Hetero":   dict(cx=12.9, cy=10.3, w=2.4, h=1.4, kind="hetero", text="Hetero-\ntrophs"),
    "Microzoo": dict(cx=15.9, cy=10.3, w=2.7, h=1.4, kind="zoo", text="Microzoopl."),
    "Mesozoo":  dict(cx=19.0, cy=10.3, w=2.7, h=1.4, kind="zoo", text="Mesozoopl."),

    # ---- PELAGIC organic pools ----
    "DOM_p":  dict(cx=23.2, cy=6.0,  w=2.6, h=1.5, kind="pool", text="DOM"),
    "Part_p": dict(cx=23.2, cy=9.3,  w=2.6, h=1.5, kind="pool", text="Particulates"),

    # ---- BENTHIC carbonate system (mirror; no air-sea pCO2) ----
    "dic_b":     dict(cx=4.4, cy=17.2, w=1.5, h=1.2, kind="gas", text="DIC"),
    "o2_b":      dict(cx=4.4, cy=18.8, w=1.5, h=1.2, kind="gas", text="O₂"),
    "dicspec_b": dict(cx=2.2, cy=20.4, w=2.2, h=2.6, kind="spec",
                      text="H₂CO₃\n⇅\nHCO₃⁻\n⇅\nCO₃²⁻"),
    "ph_b":      dict(cx=4.4, cy=20.4, w=1.6, h=1.3, kind="gas", text="pH\nΩ"),
    "ta_b":      dict(cx=2.2, cy=23.5, w=1.6, h=1.3, kind="gas", text="TA"),

    # ---- BENTHIC nutrients (Si, PO4, NO3, NH4; no Fe) ----
    "Si_b": dict(cx=6.4, cy=18.6, w=1.6, h=1.15, kind="capsule", text="Si"),
    "PO4b": dict(cx=6.4, cy=19.9, w=1.6, h=1.15, kind="capsule", text="PO₄"),
    "NO3b": dict(cx=6.4, cy=21.2, w=1.6, h=1.15, kind="capsule", text="NO₃"),
    "NH4b": dict(cx=6.4, cy=22.5, w=1.6, h=1.15, kind="capsule", text="NH₄"),

    # ---- BENTHIC zoobenthos ----
    "Suspension": dict(cx=12.9, cy=18.1, w=2.5, h=1.4, kind="zoo", text="Suspension\nFeeders"),
    "Meiobenthos": dict(cx=12.9, cy=20.2, w=2.5, h=1.4, kind="zoo", text="Meio-\nbenthos"),
    "Deposit":    dict(cx=15.8, cy=19.4, w=2.5, h=1.4, kind="zoo", text="Deposit\nFeeders"),

    # ---- BENTHIC microbes ----
    "Aerobic":   dict(cx=13.1, cy=22.1, w=2.5, h=1.2, kind="bact", text="Aerobic\nBacteria"),
    "Anaerobic": dict(cx=13.1, cy=23.7, w=2.5, h=1.2, kind="bact", text="Anaerobic\nBacteria"),

    # ---- BENTHIC organic pools ----
    "Part_b": dict(cx=23.2, cy=19.5, w=2.6, h=1.5, kind="pool", text="Particulates"),
    "DOM_b":  dict(cx=23.2, cy=24.5, w=2.6, h=1.5, kind="pool", text="DOM"),
}

# ---------------------------------------------------------------------------
# EDGES  (src_port, dst_port let the builder pick clean attachment sides)
#   route: "h" horizontal-first elbow, "v" vertical-first elbow, "s" straight
#   double: True draws arrowheads on both ends
# ---------------------------------------------------------------------------
E = lambda s, d, c, **k: dict(src=s, dst=d, color=c, **k)
EDGES = [
    # ===== PELAGIC carbonate system internal (blue equilibria) =====
    E("pco2_atm_p", "pco2_p", "blue", double=True, route="s"),          # air-sea CO2
    E("pco2_p", "dicspec_p", "blue", double=True, route="v"),           # CO2<->DIC
    E("dicspec_p", "ph_p", "blue", double=True, route="s"),             # DIC<->pH/Omega
    E("ta_p", "dicspec_p", "blue", double=True, route="v"),             # TA<->carbonate
    E("dic_p", "dicspec_p", "blue", double=True, route="v"),
    # DIC uptake by phytoplankton (photosynthesis) - two blue arrows right
    E("dic_p", "Picoph", "blue", route="h", src_port="r"),
    E("o2_p", "Picoph", "blue", route="h", src_port="r"),

    # ===== PELAGIC nutrient exchange (red) =====
    E("Fe", "Picoph", "red", route="h"),                                # nutrient uptake
    E("NO3", "Picoph", "red", route="h"),
    E("PO4", "Picoph", "red", route="h"),
    E("Bacteria_p", "PO4", "red", route="h"),                           # remineralisation
    E("Bacteria_p", "NH4", "red", route="h"),

    # ===== PELAGIC organic flows (green) =====
    E("Diatoms", "DOM_p", "green", route="v"),                          # exudation/lysis
    E("Diatoms", "Part_p", "green", route="v"),                         # mortality -> POM
    E("DOM_p", "Bacteria_p", "green", route="h"),                       # DOM uptake
    E("DOM_p", "Part_p", "green", route="v"),                           # aggregation

    # ===== PELAGIC grazing (black) =====
    E("Picoph", "Hetero", "black", route="v"),                          # microbial loop graze
    E("Microph", "Microzoo", "black", route="v"),
    E("Diatoms", "Mesozoo", "black", route="v"),
    E("Bacteria_p", "Hetero", "black", route="v"),
    E("Hetero", "Microzoo", "black", route="s"),                        # trophic transfer
    E("Microzoo", "Mesozoo", "black", route="s"),
    E("Mesozoo", "Part_p", "green", route="h"),                         # egestion/mortality

    # ===== Sinking flux (grey, thick) =====
    E("Part_p", "Part_b", "grey", route="s", thick=True),

    # ===== Pelagic<->Benthic solute exchange across the sediment-water =====
    # interface (short blue dashed double arrows; original draws "curl" symbols)
    dict(raw=True, x1=3.0, y1=15.2, x2=3.0, y2=16.8, color="blue",
         double=True, dash=True),
    dict(raw=True, x1=4.4, y1=15.2, x2=4.4, y2=16.8, color="blue",
         double=True, dash=True),

    # ===== BENTHIC carbonate internal (blue) =====
    E("dic_b", "dicspec_b", "blue", double=True, route="v"),
    E("dicspec_b", "ph_b", "blue", double=True, route="s"),
    E("ta_b", "dicspec_b", "blue", double=True, route="v"),

    # ===== BENTHIC nutrient exchange (red/blue) - aerobic vs anaerobic distinct =====
    E("Aerobic", "Si_b", "blue", route="h"),                           # silicate dissolution
    E("Aerobic", "PO4b", "red", route="h", double=True),               # P exchange (two-way)
    E("Anaerobic", "NO3b", "red", route="h"),                          # denitrification
    E("Anaerobic", "NH4b", "red", route="h"),                          # ammonification

    # ===== BENTHIC organic flows (green) =====
    E("Part_b", "DOM_b", "green", route="v"),
    E("Part_b", "Aerobic", "green", route="h"),                        # POM -> bacteria
    E("DOM_b", "Anaerobic", "green", route="h"),                       # DOM -> bacteria
    E("DOM_b", "Part_b", "green", route="v"),

    # ===== BENTHIC grazing (black) =====
    E("Part_b", "Suspension", "black", route="h"),
    E("Part_b", "Deposit", "black", route="h"),
    E("Meiobenthos", "Deposit", "black", route="s"),

    # ===== Sinking from pelagic consumers/particulates already via Part_p =====
]
