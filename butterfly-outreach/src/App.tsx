import { useState, useEffect, useRef } from "react";
import {
  ChevronDown, ChevronUp, ExternalLink, Mail, BookOpen,
  Globe, ArrowDown, Layers, BarChart2, GitBranch, FlaskConical,
} from "lucide-react";
import { figures } from "./assets/figures";

// ─── Inject Google Fonts at runtime (keeps html-inline from choking) ─────────
function injectFonts() {
  const id = "gf-butterfly";
  if (document.getElementById(id)) return;
  const urls = [
    "https://fonts.googleapis.com",
    "https://fonts.gstatic.com",
  ];
  urls.forEach((href, i) => {
    const l = document.createElement("link");
    l.rel = "preconnect";
    l.href = href;
    if (i === 1) l.crossOrigin = "anonymous";
    document.head.appendChild(l);
  });
  const l = document.createElement("link");
  l.id = id;
  l.rel = "stylesheet";
  l.href =
    "https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&family=Inter:wght@300;400;500;600&display=swap";
  document.head.appendChild(l);
}

// ─── Design tokens ────────────────────────────────────────────────────────────
const T = {
  // backgrounds
  bg0:    "#ffffff",           // pure white (hero, main)
  bg1:    "#f2fbf8",           // barely-emerald section tint
  bg2:    "#e6f7f2",           // light card surface
  bgDark: "#053d22",           // deep forest (contrast section)
  bgDark2:"#074d2b",           // slightly lifted forest

  // borders
  border:     "rgba(0,168,107,0.15)",
  borderDark: "rgba(255,255,255,0.10)",

  // primary spectrum — emerald → teal → sky
  emerald:  "#00a86b",   // radiant emerald green
  emeraldL: "#00c48c",   // lighter emerald
  teal:     "#00b8a9",   // luminous teal
  tealL:    "#00d4c8",   // glowing teal
  sky:      "#0ea5e9",   // sky blue
  skyL:     "#38bdf8",   // light sky
  cyan:     "#22d3ee",   // soft glowing cyan

  // accent
  coral:   "#ff6b6b",   // warm coral — key findings & CTAs

  // text on light bg
  text0:  "#052e1c",    // deep forest — primary
  text1:  "#2d6b50",    // secondary
  text2:  "#6aaa8f",    // muted

  // text on dark bg
  dtxt0:  "#e8fff6",
  dtxt1:  "#8dd4b5",
  dtxt2:  "#3d8f66",

  mono:    "'JetBrains Mono', 'Courier New', monospace",
  display: "'Outfit', system-ui, sans-serif",
  body:    "'Inter', system-ui, sans-serif",
};

// ─── Dot-grid SVG backgrounds ─────────────────────────────────────────────────
// Light sections: emerald dots
const DOT_GRID = `url("data:image/svg+xml,%3Csvg width='28' height='28' viewBox='0 0 28 28' xmlns='http://www.w3.org/2000/svg'%3E%3Ccircle cx='1' cy='1' r='1' fill='%2300a86b' fill-opacity='0.12'/%3E%3C/svg%3E")`;
// Dark sections: white dots
const DOT_GRID_DARK = `url("data:image/svg+xml,%3Csvg width='28' height='28' viewBox='0 0 28 28' xmlns='http://www.w3.org/2000/svg'%3E%3Ccircle cx='1' cy='1' r='1' fill='%23ffffff' fill-opacity='0.06'/%3E%3C/svg%3E")`;

// ─── Scroll hook ──────────────────────────────────────────────────────────────
function useScrolled(threshold = 60) {
  const [v, setV] = useState(false);
  useEffect(() => {
    const fn = () => setV(window.scrollY > threshold);
    window.addEventListener("scroll", fn, { passive: true });
    return () => window.removeEventListener("scroll", fn);
  }, [threshold]);
  return v;
}

// ─── Fade-in hook ─────────────────────────────────────────────────────────────
function useFadeIn() {
  const ref = useRef<HTMLDivElement>(null);
  const [visible, setVisible] = useState(false);
  useEffect(() => {
    const obs = new IntersectionObserver(
      ([e]) => { if (e.isIntersecting) { setVisible(true); obs.disconnect(); } },
      { threshold: 0.08 }
    );
    if (ref.current) obs.observe(ref.current);
    return () => obs.disconnect();
  }, []);
  return { ref, visible };
}

// ─── Section ──────────────────────────────────────────────────────────────────
function Section({ id, className = "", style = {}, children }: {
  id?: string; className?: string; style?: React.CSSProperties; children: React.ReactNode;
}) {
  const { ref, visible } = useFadeIn();
  return (
    <section
      id={id} ref={ref} style={style}
      className={`transition-all duration-700 ease-out ${visible ? "opacity-100 translate-y-0" : "opacity-0 translate-y-6"} ${className}`}
    >
      {children}
    </section>
  );
}

// ─── Label chip ───────────────────────────────────────────────────────────────
function Label({ children }: { children: React.ReactNode }) {
  return (
    <p style={{ fontFamily: T.mono, color: T.teal, fontSize: 11, letterSpacing: "0.15em" }}
       className="uppercase mb-3 font-medium">
      {children}
    </p>
  );
}

// ─── Section heading ──────────────────────────────────────────────────────────
function Heading({ children }: { children: React.ReactNode }) {
  return (
    <h2 style={{ fontFamily: T.display, color: T.text0, letterSpacing: "-0.03em" }}
        className="text-3xl md:text-4xl font-bold leading-tight mb-2">
      {children}
    </h2>
  );
}

// ─── Monospace data tag ───────────────────────────────────────────────────────
function DataTag({ children, accent = T.teal }: { children: React.ReactNode; accent?: string }) {
  return (
    <code style={{
      fontFamily: T.mono, fontSize: 12, color: accent,
      background: `${accent}14`, border: `1px solid ${accent}30`,
      padding: "2px 8px", borderRadius: 3,
    }}>
      {children}
    </code>
  );
}

// ─── Stat block ───────────────────────────────────────────────────────────────
function Stat({ value, label, accent = T.teal }: { value: string; label: string; accent?: string }) {
  return (
    <div className="flex flex-col gap-1">
      <span style={{ fontFamily: T.mono, color: accent, fontSize: 32, fontWeight: 600, lineHeight: 1 }}>
        {value}
      </span>
      <span style={{ fontFamily: T.body, color: T.text2, fontSize: 11, letterSpacing: "0.1em" }}
            className="uppercase">
        {label}
      </span>
    </div>
  );
}

// ─── Accordion ────────────────────────────────────────────────────────────────
// Parses "Panel a — Description text" into a badge + label pair.
function Accordion({ title, children }: { title: string; children: React.ReactNode }) {
  const [open, setOpen] = useState(false);
  const m = title.match(/^(Panel [a-z\d]+)\s*[—–-]+\s*(.+)$/i);
  const badge = m ? m[1] : null;
  const label = m ? m[2] : title;

  return (
    <div style={{ borderTop: `1px solid ${T.border}` }}>
      <button
        onClick={() => setOpen(!open)}
        className="w-full flex items-center justify-between gap-4 py-3.5 text-left"
        style={{ background: "transparent" }}
      >
        <div className="flex items-center gap-3 min-w-0">
          {badge && (
            <span style={{
              fontFamily: T.mono, fontSize: 10, fontWeight: 600, letterSpacing: "0.1em",
              color: T.emerald, background: `${T.emerald}14`,
              border: `1px solid ${T.emerald}35`,
              padding: "2px 9px", flexShrink: 0, whiteSpace: "nowrap",
            }}>
              {badge.toUpperCase()}
            </span>
          )}
          <span style={{
            fontFamily: T.body, fontWeight: 500, fontSize: 13.5, lineHeight: 1.45,
            color: open ? T.emerald : T.text0,
            transition: "color 0.18s",
          }}>
            {label}
          </span>
        </div>
        <span style={{ color: T.emerald, flexShrink: 0 }}>
          {open ? <ChevronUp size={14} /> : <ChevronDown size={14} />}
        </span>
      </button>

      <div className="overflow-hidden transition-all duration-300"
           style={{ maxHeight: open ? 1000 : 0 }}>
        <div style={{
          fontFamily: T.body, color: T.text0,
          fontSize: 14, lineHeight: 1.82,
          paddingBottom: 22,
          paddingLeft: 16,
          borderLeft: `2px solid ${T.emerald}25`,
          marginLeft: 1,
          marginBottom: 4,
        }}>
          {children}
        </div>
      </div>
    </div>
  );
}

// ─── Finding card ─────────────────────────────────────────────────────────────
function FindingCard({ num, title, bullets, accent, tags }: {
  num: string; title: string; bullets: string[]; accent: string; tags: string[];
}) {
  return (
    <div
      className="relative p-7 flex flex-col gap-5"
      style={{
        background: T.bg2,
        border: `1px solid ${T.border}`,
        borderTop: `2px solid ${accent}`,
      }}
    >
      <div className="flex items-start justify-between gap-4">
        <span style={{ fontFamily: T.mono, color: T.text2, fontSize: 11 }}>F-{num}</span>
        <div className="flex gap-2 flex-wrap justify-end">
          {tags.map(t => <DataTag key={t} accent={accent}>{t}</DataTag>)}
        </div>
      </div>
      <h3 style={{ fontFamily: T.display, color: T.text0, fontSize: 17, fontWeight: 600, lineHeight: 1.4 }}>
        {title}
      </h3>
      <ul className="space-y-2.5">
        {bullets.map((b, i) => (
          <li key={i} className="flex gap-3" style={{ fontFamily: T.body, color: T.text1, fontSize: 14, lineHeight: 1.65 }}>
            <span style={{ color: accent, fontFamily: T.mono, fontSize: 12, marginTop: 2, flexShrink: 0 }}>—</span>
            {b}
          </li>
        ))}
      </ul>
    </div>
  );
}

// ─── Method step ──────────────────────────────────────────────────────────────
function Step({ n, icon, title, body }: { n: number; icon: React.ReactNode; title: string; body: string }) {
  return (
    <div className="flex gap-5">
      <div className="flex flex-col items-center gap-0">
        <div
          style={{
            width: 36, height: 36, border: `1px solid ${n % 2 ? T.teal : T.coral}`,
            color: n % 2 ? T.teal : T.coral, fontFamily: T.mono, fontSize: 13, flexShrink: 0,
          }}
          className="flex items-center justify-center font-semibold"
        >
          {String(n).padStart(2, "0")}
        </div>
        {n < 5 && <div style={{ width: 1, flex: 1, background: T.border, minHeight: 40 }} className="my-1" />}
      </div>
      <div className="pb-10">
        <div style={{ color: T.text2, marginBottom: 6 }}>{icon}</div>
        <h4 style={{ fontFamily: T.display, color: T.text0, fontSize: 15, fontWeight: 600, marginBottom: 6 }}>
          {title}
        </h4>
        <p style={{ fontFamily: T.body, color: T.text1, fontSize: 13.5, lineHeight: 1.7 }}>{body}</p>
      </div>
    </div>
  );
}

// ─── Figure container ─────────────────────────────────────────────────────────
function FigureCard({ num, accent, title, sub, src, alt, children }: {
  num: string; accent: string; title: string; sub: string;
  src: string; alt: string; children: React.ReactNode;
}) {
  return (
    <div style={{ background: T.bg1, border: `1px solid ${T.border}` }}>
      {/* Header bar */}
      <div style={{ background: T.bg2, borderBottom: `1px solid ${T.border}`, padding: "14px 24px" }}
           className="flex items-center justify-between">
        <div className="flex items-center gap-4">
          <span style={{
            fontFamily: T.mono, fontSize: 11, color: accent,
            border: `1px solid ${accent}40`, padding: "2px 10px",
          }}>
            FIG.{num}
          </span>
          <span style={{ fontFamily: T.display, color: T.text0, fontWeight: 600, fontSize: 15 }}>{title}</span>
        </div>
        <span style={{ fontFamily: T.mono, color: T.text2, fontSize: 11 }}>{sub}</span>
      </div>
      {/* Figure image */}
      <div style={{ background: "#fff", padding: "0" }}>
        <img src={src} alt={alt} className="w-full h-auto block" style={{ display: "block" }} />
      </div>
      {/* Accordions */}
      <div style={{
        padding: "4px 24px 16px",
        borderTop: `1px solid ${T.border}`,
        background: T.bg0,
      }}>
        {children}
      </div>
    </div>
  );
}

// ─── Nav ──────────────────────────────────────────────────────────────────────
const NAV = [
  { l: "Findings", h: "#findings" },
  { l: "Figures", h: "#figures" },
  { l: "Methods", h: "#methods" },
  { l: "Conservation", h: "#conservation" },
  { l: "About", h: "#about" },
];

function Nav() {
  const scrolled = useScrolled();
  return (
    <nav style={{
      position: "fixed", top: 0, left: 0, right: 0, zIndex: 50,
      background: scrolled ? "rgba(255,255,255,0.95)" : "transparent",
      backdropFilter: scrolled ? "blur(16px) saturate(180%)" : "none",
      borderBottom: scrolled ? `1px solid ${T.border}` : "none",
      boxShadow: scrolled ? "0 1px 24px rgba(0,168,107,0.08)" : "none",
      transition: "all 0.3s",
    }}>
      <div className="max-w-6xl mx-auto px-6 h-14 flex items-center justify-between">
        <span style={{ fontFamily: T.mono, color: T.emerald, fontSize: 12, letterSpacing: "0.08em" }}>
          LEPIDOPTERA / MIGRATION / GLOBAL
        </span>
        <div className="hidden md:flex gap-8">
          {NAV.map(n => (
            <a key={n.h} href={n.h} style={{ fontFamily: T.body, color: T.text1, fontSize: 13 }}
               className="transition-colors hover:opacity-70">
              {n.l}
            </a>
          ))}
        </div>
        <a href="#about" style={{
          fontFamily: T.mono, fontSize: 11, color: "#fff",
          background: T.emerald, padding: "6px 16px", letterSpacing: "0.08em",
        }}>
          READ PAPER
        </a>
      </div>
    </nav>
  );
}

// ─── App ──────────────────────────────────────────────────────────────────────
export default function App() {
  useEffect(() => { injectFonts(); }, []);
  return (
    <div style={{ background: T.bg0, color: T.text0, fontFamily: T.body, minHeight: "100vh" }}>
      <Nav />

      {/* ── HERO ──────────────────────────────────────────────────────── */}
      <header style={{
        minHeight: "100vh", display: "flex", alignItems: "center", justifyContent: "center",
        background: "linear-gradient(150deg, #ffffff 0%, #f0fbf7 45%, #e2f8f2 100%)",
        backgroundImage: `${DOT_GRID}, linear-gradient(150deg, #ffffff 0%, #f0fbf7 45%, #e2f8f2 100%)`,
        backgroundSize: "28px 28px, cover",
        position: "relative", overflow: "hidden",
        padding: "80px 24px 60px",
      }}>
        {/* Radial glow — bioluminescent emerald */}
        <div style={{
          position: "absolute", inset: 0, pointerEvents: "none",
          background: "radial-gradient(ellipse 70% 55% at 55% 38%, rgba(0,196,140,0.18) 0%, rgba(14,165,233,0.06) 55%, transparent 80%)",
        }} />
        {/* Corner coordinate labels */}
        <span style={{ position: "absolute", top: 80, left: 24, fontFamily: T.mono, color: T.text2, fontSize: 11 }}>
          90°N ↓ 90°S
        </span>
        <span style={{ position: "absolute", top: 80, right: 24, fontFamily: T.mono, color: T.text2, fontSize: 11 }}>
          180°W ↔ 180°E
        </span>

        <div style={{ maxWidth: 820, width: "100%", zIndex: 1 }}>
          {/* Tag row */}
          <div className="flex gap-3 flex-wrap mb-8">
            {["n = 426 species", "4 seasons", "5 families", "Global scale"].map(t => (
              <span key={t} style={{
                fontFamily: T.mono, fontSize: 11, color: T.teal, letterSpacing: "0.1em",
                border: `1px solid ${T.teal}30`, padding: "3px 12px",
              }}>
                {t.toUpperCase()}
              </span>
            ))}
          </div>

          {/* Main title */}
          <h1 style={{
            fontFamily: T.display, color: T.text0,
            fontSize: "clamp(36px, 6vw, 72px)",
            fontWeight: 700, lineHeight: 1.08, letterSpacing: "-0.04em",
            marginBottom: 28,
          }}>
            Why Do Migratory Butterflies<br />
            <span style={{ color: T.teal }}>Follow Their Ancestors?</span>
          </h1>

          <p style={{ fontFamily: T.body, color: T.text1, fontSize: 18, lineHeight: 1.7, maxWidth: 620, marginBottom: 48 }}>
            A global macroecological analysis of 426 species reveals that evolutionary heritage — not just climate — structures where butterflies migrate across seasons and continents.
          </p>

          {/* Stats row */}
          <div style={{
            display: "grid", gridTemplateColumns: "repeat(4, 1fr)",
            gap: 0,
            border: `1px solid ${T.border}`,
            marginBottom: 40, maxWidth: 680,
          }}>
            {[
              { v: "89%",  l: "Variance from phylogeny",  a: T.teal },
              { v: "12°N", l: "Global richness peak",     a: T.coral },
              { v: "0.727",l: "Rapoport β coefficient",   a: T.sky },
              { v: "3.311",l: "Grassland impact β",       a: T.cyan },
            ].map((s, i) => (
              <div key={s.l} style={{
                padding: "20px 20px",
                borderRight: i < 3 ? `1px solid ${T.border}` : "none",
                background: T.bg1,
              }}>
                <Stat value={s.v} label={s.l} accent={s.a} />
              </div>
            ))}
          </div>

          <div className="flex gap-4 flex-wrap">
            <a href="#findings" style={{
              fontFamily: T.mono, fontSize: 12, letterSpacing: "0.1em",
              color: "#fff", background: T.emerald, padding: "12px 28px",
              display: "flex", alignItems: "center", gap: 8,
            }}>
              EXPLORE FINDINGS <ArrowDown size={14} />
            </a>
            <a href="#figures" style={{
              fontFamily: T.mono, fontSize: 12, letterSpacing: "0.1em",
              color: T.emerald, border: `1.5px solid ${T.emerald}`, padding: "12px 28px",
            }}>
              VIEW FIGURES
            </a>
          </div>
        </div>

        {/* Bottom scroll cue */}
        <div style={{ position: "absolute", bottom: 28, left: "50%", transform: "translateX(-50%)", color: T.text2 }}
             className="animate-bounce">
          <ChevronDown size={20} />
        </div>
      </header>

      {/* ── FINDINGS ──────────────────────────────────────────────────── */}
      <div style={{ background: T.bg1, backgroundImage: DOT_GRID, backgroundSize: "28px 28px" }}>
        <Section id="findings" className="max-w-6xl mx-auto px-6 py-24">
          <Label>Core Discoveries</Label>
          <Heading>Four Results That Reframe<br />Butterfly Migration</Heading>
          <p style={{ fontFamily: T.body, color: T.text1, fontSize: 15, maxWidth: 560, marginBottom: 48, marginTop: 12, lineHeight: 1.7 }}>
            Integrating species distribution models, morphological traits, and Bayesian phylogenetic models across all major migratory families.
          </p>

          <div style={{ display: "grid", gridTemplateColumns: "repeat(2, 1fr)", gap: 1, background: T.border }}>
            <FindingCard
              num="01" accent={T.teal}
              title="Tropical richness peaks at 12°N — not in the temperate zone"
              tags={["β richness", "12°N peak"]}
              bullets={[
                "Mean migratory richness peaks at 12°N (16.95 species per cell), contradicting the expectation that migration is a temperate-zone phenomenon.",
                "The Northern Hemisphere peak exceeds the Southern Hemisphere maximum (7.73 species near 30°S) by more than 2×.",
                "In spring (S2→S1), 55.2% of the global migratory range records a net species gain as butterflies advance poleward.",
              ]}
            />
            <FindingCard
              num="02" accent={T.coral}
              title="Rapoport's Rule confirmed — tropics contract seasonal ranges"
              tags={["β = −0.727", "p < 0.001"]}
              bullets={[
                "Tropical-restricted species have significantly smaller seasonal ranges: every 10% increase in tropical occupancy corresponds to a 0.73-unit drop in log₁₀ range size.",
                "Temperate and boreal species maintain ranges roughly 5× larger than their tropical counterparts.",
                "Seasonal range size peaks during Northern Hemisphere summer as species colonise high-latitude zones.",
              ]}
            />
            <FindingCard
              num="03" accent={T.sky}
              title="Wing size is a latitudinal lever — but only above 18°"
              tags={["β lat×WS = 0.062", "p = 0.008"]}
              bullets={[
                "Wingspan decreases as absolute latitude decreases (β = −0.002, p < 0.001) — larger-winged species occupy higher latitudes.",
                "At 18°N/S the slope between wingspan and range size converges to zero; above that inflection it becomes strongly positive.",
                "Morphological dispersal capacity only translates to range expansion in seasonally severe environments.",
              ]}
            />
            <FindingCard
              num="04" accent={T.cyan}
              title="Evolutionary ancestry accounts for 89% of range variance"
              tags={["H² = 0.89", "λ ≈ 0.90"]}
              bullets={[
                "Bayesian phylogenetic models (BPMM) show H² = 0.89 — closely related species share range size regardless of contemporary climate.",
                "Grassland conversion is the single strongest anthropogenic predictor of hotspot loss (OLS β = −3.311, p < 0.001).",
                "High dispersal capacity preserves ancestral niches rather than enabling in-situ climate adaptation.",
              ]}
            />
          </div>
        </Section>
      </div>

      {/* ── FIGURES ───────────────────────────────────────────────────── */}
      <div style={{ background: T.bg0 }}>
        <Section id="figures" className="max-w-6xl mx-auto px-6 py-24">
          <Label>Data Visualisations</Label>
          <Heading>The Empirical Evidence</Heading>
          <p style={{ fontFamily: T.body, color: T.text1, fontSize: 15, maxWidth: 540, marginBottom: 48, marginTop: 12, lineHeight: 1.7 }}>
            Three figures from the manuscript with annotated statistical interpretations.
          </p>

          <div className="flex flex-col gap-6">

            <FigureCard
              num="01" accent={T.teal}
              title="Global Seasonal Migration Patterns"
              sub="Richness · Latitudinal profile · Seasonal flux"
              src={figures.Figure1}
              alt="Figure 1 – Global richness and seasonal dynamics"
            >
              <Accordion title="Panel a — How to read the world map">
                <p>Each 0.5° grid cell is coloured by the count of migratory species present as seasonal switchers. Deep teal clusters in Central America, sub-Saharan Africa, and South/Southeast Asia mark global hotspots. Notably, the Amazon and Congo basins appear muted — not from species absence, but because many tropical residents remain present across all seasons and are excluded from the migratory subset by design.</p>
              </Accordion>
              <Accordion title="Panel b — The 12°N anomaly">
                <p>The latitudinal profile shows mean migratory richness (y) against latitude (x). The spike at <DataTag>12°N = 16.95 spp</DataTag> is the study's first major result: subtropical zones — not high-latitude temperate belts — host the densest migrant assemblages. A secondary peak near 30°S confirms hemispheric asymmetry, likely amplified by data gaps in southern tropical regions.</p>
              </Accordion>
              <Accordion title="Panel c — Reading the seasonal flux maps">
                <p>Blue pixels = net species arrivals; red = departures. The <DataTag>S2→S1</DataTag> (spring) panel shows a poleward surge: 55.2% of the migratory range gains species. The <DataTag>S4→S3</DataTag> (autumn) panel reverses this, with 55.7% of the range recording net losses as individuals retreat equatorward to overwintering habitats.</p>
              </Accordion>
            </FigureCard>

            <FigureCard
              num="02" accent={T.coral}
              title="Morphological, Geographic & Environmental Drivers"
              sub="Range size · Wingspan · Interaction · Forest plot"
              src={figures.Figure2}
              alt="Figure 2 – Morphology and environmental predictors"
            >
              <Accordion title="Panel a — Rapoport's Rule quantified">
                <p>Scatter shows log₁₀ range size (y) vs. mean tropical occupancy (x), coloured by temperature seasonality <DataTag>Bio4</DataTag>. The LMM slope of <DataTag>β = −0.727</DataTag> confirms Rapoport's Rule: a species spending its full range in the tropics has a range ~5× smaller than one centred in temperate latitudes. The colour gradient adds a third dimension: warm-coloured (high seasonality) points cluster among large-range species, linking climatically variable environments to broad distributions.</p>
              </Accordion>
              <Accordion title="Panel c — The wingspan × latitude interaction (key result)">
                <p>This interaction plot is the study's morphological centrepiece. The slope between wingspan and range size is plotted against absolute latitude. At <DataTag>18°N/S</DataTag> the slope is ~0. Below 18° (tropics), wingspan is irrelevant to range size. Above 18°, a larger wingspan increasingly predicts a larger range — suggesting that dispersal morphology only becomes a binding constraint where seasonal climate creates genuine barriers to range expansion.</p>
              </Accordion>
              <Accordion title="Panel d — Forest plot: what moves the needle on hotspot richness?">
                <ul style={{ listStyle: "none", padding: 0 }} className="space-y-2">
                  <li><DataTag accent={T.coral}>Grassland conversion β = −3.311</DataTag> — strongest negative predictor. Converting native grasslands eliminates essential nectar sources and larval habitats across migratory corridors.</li>
                  <li className="mt-2"><DataTag accent={T.sky}>Shrubland β = −3.117</DataTag> and <DataTag accent={T.sky}>Forest β = −2.408</DataTag> — similarly negative, reinforcing that natural land cover is structurally critical.</li>
                  <li className="mt-2"><DataTag accent={T.teal}>Temperature seasonality β = −0.689</DataTag> — counter-intuitively negative: extreme seasonality is too harsh even for mobile species.</li>
                  <li className="mt-2"><DataTag>HII β = +0.137</DataTag> — slight positive, likely a sampling artefact from citizen-science density near human populations.</li>
                </ul>
              </Accordion>
            </FigureCard>

            <FigureCard
              num="03" accent={T.cyan}
              title="Phylogenetic Signal in Range Size"
              sub="Circular phylogeny · 247 species · 5 families"
              src={figures.Figure3}
              alt="Figure 3 – Phylogenetic framework and range size distribution"
            >
              <Accordion title="How to read the circular phylogeny">
                <p>Each leaf tip represents a species; branches cluster into five families: Pieridae, Hesperiidae, Papilionidae, Nymphalidae, and Lycaenidae. Branch colours encode log-transformed range size — warm hues (yellow → red) = large ranges; cool hues (blue → violet) = small ranges. The inset histogram shows the overall distribution of range sizes across the sample.</p>
              </Accordion>
              <Accordion title="What H² = 0.89 actually means">
                <p>If range size were randomly distributed across the phylogeny, colours would be mixed throughout the tree. Instead, colour patches cluster within clades — closely related species share similar range sizes. The BPMM quantifies this as <DataTag accent={T.cyan}>H² = 0.89</DataTag> (95% CI [0.86, 0.91]): after accounting for climate, season, and morphology, 89% of the remaining variance in range size is phylogenetically structured. Evolutionary lineage is the dominant predictor of contemporary migratory range — more than any measured environmental variable.</p>
              </Accordion>
              <Accordion title="Phylogenetic niche conservatism vs. constraint">
                <p>A common misreading: H² = 0.89 does not mean these butterflies are passive victims of their ancestry. The study argues that high dispersal capacity and phylogenetic niche conservatism (PNC) are mutually reinforcing. Mobility allows butterflies to <em>follow</em> their ancestral climatic niches across space, removing the selective pressure to adapt in place. The evolutionary signal reflects strategic spatial tracking, not evolutionary stasis.</p>
              </Accordion>
            </FigureCard>

          </div>
        </Section>
      </div>

      {/* ── METHODS ───────────────────────────────────────────────────── */}
      <div style={{ background: T.bg1, borderTop: `1px solid ${T.border}` }}>
        <Section id="methods" className="max-w-6xl mx-auto px-6 py-24">
          <div className="grid md:grid-cols-2 gap-16 items-start">
            <div>
              <Label>Analytical Pipeline</Label>
              <Heading>Multi-Method,<br />Multi-Scale</Heading>
              <p style={{ fontFamily: T.body, color: T.text1, fontSize: 15, lineHeight: 1.7, marginTop: 16, marginBottom: 32 }}>
                A hierarchical inference framework integrating species distribution modelling, spatial statistics, and Bayesian evolutionary models across four seasons and five families.
              </p>
              <div style={{
                background: T.bg2, border: `1px solid ${T.border}`, padding: "16px 20px",
                display: "flex", alignItems: "center", gap: 10,
              }}>
                <FlaskConical size={15} style={{ color: T.teal, flexShrink: 0 }} />
                <span style={{ fontFamily: T.mono, color: T.text1, fontSize: 12 }}>
                  R 4.4.2 · Python 3.13 · brms · lme4 · statsmodels · terra
                </span>
              </div>
            </div>
            <div>
              <Step n={1} icon={<Globe size={14} />}
                title="Seasonal Suitability Maps"
                body="Habitat suitability maps for 426 species (Chowdhury et al. 2025), four seasons, all AUC > 0.70. Seasonal 'switchers' identified: grid cells where species presence < global seasonal maximum." />
              <Step n={2} icon={<BarChart2 size={14} />}
                title="Trait & Environmental Metrics"
                body="Range size (km²), tropical occupancy from suitability maps; wingspan from LepTraits 1.1.0 (377 spp.); WorldClim Bio4/Bio15; CLGS land cover; Human Influence Index." />
              <Step n={3} icon={<Layers size={14} />}
                title="Rapoport & Bergmann Tests (LMM)"
                body="Linear Mixed Models (lme4): range size ~ tropical occupancy + season, species as random effect; interaction model for wingspan × |latitude| to test dispersal-latitude dependency." />
              <Step n={4} icon={<BarChart2 size={14} />}
                title="Hotspot Drivers (OLS + GLM)"
                body="OLS with 100-iteration stratified resampling; GLM Poisson with spatial thinning at 50/100/200 km (Moran's I reduced from 0.485 → ~0). All predictors z-score standardised." />
              <Step n={5} icon={<GitBranch size={14} />}
                title="Bayesian Phylogenetic Models (BPMM)"
                body="brms + Kawahara et al. 2023 phylogeny; 247 matched species; 4 HMC chains × 6,000 iterations (2,000 warmup); R̂ < 1.01. Phylogenetic heritability H² partitions evolutionary vs. environmental variance." />
            </div>
          </div>
        </Section>
      </div>

      {/* ── CONSERVATION ─────────────────────────────────────────────── */}
      <div style={{
        background: T.bgDark,
        backgroundImage: DOT_GRID_DARK,
        backgroundSize: "28px 28px",
      }}>
        <Section id="conservation" className="max-w-6xl mx-auto px-6 py-24">
          {/* label override for dark bg */}
          <p style={{ fontFamily: T.mono, color: T.emeraldL, fontSize: 11, letterSpacing: "0.15em", marginBottom: 12 }}
             className="uppercase font-medium">
            Conservation Implications
          </p>
          <h2 style={{ fontFamily: T.display, color: T.dtxt0, letterSpacing: "-0.03em", fontSize: "clamp(28px,4vw,42px)", fontWeight: 700, marginBottom: 8 }}>
            High Mobility ≠ High Resilience
          </h2>
          <p style={{ fontFamily: T.body, color: T.dtxt1, fontSize: 15, maxWidth: 580, lineHeight: 1.7, margin: "16px 0 48px" }}>
            Migratory butterflies use dispersal to maintain ancestral climatic niches — but this strategy assumes those niches persist.
          </p>

          <div style={{
            display: "grid", gridTemplateColumns: "repeat(3, 1fr)",
            border: `1px solid ${T.borderDark}`, gap: 0, background: T.borderDark,
          }}>
            {[
              {
                n: "01", accent: T.coral, icon: "◈",
                title: "Protect Grassland Corridors",
                body: "Grassland conversion carries the highest single negative effect on migratory richness (β = −3.311). Native grasslands along flyways provide irreplaceable nectar and larval habitat — their loss is not compensated by adjacent habitats.",
              },
              {
                n: "02", accent: T.skyL, icon: "◈",
                title: "Climate Tracking ≠ Adaptation",
                body: "These butterflies spatially track ancestral climatic envelopes rather than genetically adapting to novel climates. When seasonal envelopes shift beyond reachable space, dispersal capacity provides no buffer.",
              },
              {
                n: "03", accent: T.tealL, icon: "◈",
                title: "Phylogenetic Triage",
                body: "H² = 0.89 means phylogenetic position is the best predictor of range constraints. Lineages with concentrated evolutionary relatedness can be prioritised for conservation assessment before empirical range data are available.",
              },
            ].map(c => (
              <div key={c.n} style={{ background: T.bgDark2, padding: "28px 24px", borderTop: `2px solid ${c.accent}` }}>
                <div style={{ fontFamily: T.mono, color: c.accent, fontSize: 22, marginBottom: 12 }}>{c.icon}</div>
                <p style={{ fontFamily: T.mono, color: c.accent, fontSize: 11, letterSpacing: "0.1em", marginBottom: 10 }}>
                  IMPLICATION {c.n}
                </p>
                <h3 style={{ fontFamily: T.display, color: T.dtxt0, fontWeight: 600, fontSize: 16, marginBottom: 12 }}>
                  {c.title}
                </h3>
                <p style={{ fontFamily: T.body, color: T.dtxt1, fontSize: 14, lineHeight: 1.7 }}>{c.body}</p>
              </div>
            ))}
          </div>
        </Section>
      </div>

      {/* ── ABOUT ─────────────────────────────────────────────────────── */}
      <div style={{ background: T.bg1, borderTop: `1px solid ${T.border}` }}>
        <Section id="about" className="max-w-6xl mx-auto px-6 py-24">
          <div className="grid md:grid-cols-2 gap-16 items-start">
            <div>
              <Label>About the Study</Label>
              <Heading>Research &amp; Contact</Heading>
              <p style={{ fontFamily: T.body, color: T.text1, fontSize: 15, lineHeight: 1.7, margin: "16px 0 32px" }}>
                This study combines global biodiversity informatics, morphological trait databases, and Bayesian evolutionary modelling to test macroecological rules in migratory insects at an unprecedented global scale.
              </p>

              <div className="flex flex-col gap-3">
                {[
                  { label: "Full Manuscript", sub: "Preprint / peer-reviewed article", icon: <BookOpen size={15} /> },
                  { label: "Data & Code", sub: "GitHub · R + Python scripts", icon: <ExternalLink size={15} /> },
                  { label: "Collaboration", sub: "lhy09260416@gmail.com", icon: <Mail size={15} /> },
                ].map(item => (
                  <div key={item.label} style={{
                    border: `1px solid ${T.border}`, background: T.bg2,
                    padding: "14px 18px", display: "flex", alignItems: "center", gap: 14,
                    cursor: "pointer",
                    transition: "border-color 0.2s",
                  }}
                    className="hover:border-teal-400 group"
                  >
                    <div style={{ color: T.teal }}>{item.icon}</div>
                    <div style={{ flex: 1 }}>
                      <p style={{ fontFamily: T.display, color: T.text0, fontSize: 14, fontWeight: 600 }}>{item.label}</p>
                      <p style={{ fontFamily: T.mono, color: T.text2, fontSize: 11 }}>{item.sub}</p>
                    </div>
                    <ExternalLink size={13} style={{ color: T.text2 }} />
                  </div>
                ))}
              </div>
            </div>

            {/* Citation card */}
            <div style={{
              background: T.bg2, border: `1px solid ${T.border}`,
              borderTop: `2px solid ${T.teal}`,
            }}>
              <div style={{ padding: "20px 24px", borderBottom: `1px solid ${T.border}` }}>
                <span style={{ fontFamily: T.mono, color: T.teal, fontSize: 11, letterSpacing: "0.1em" }}>
                  MANUSCRIPT TITLE
                </span>
              </div>
              <div style={{ padding: "24px" }}>
                <p style={{ fontFamily: T.display, color: T.text0, fontSize: 16, fontWeight: 600, lineHeight: 1.5, marginBottom: 20 }}>
                  Seasonal climate and evolutionary history shape global hotspots of migratory butterflies
                </p>
                <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 16, marginBottom: 28 }}>
                  {[
                    { l: "Species", v: "426" },
                    { l: "Seasons", v: "4" },
                    { l: "Phylogenetic signal", v: "H² = 0.89" },
                    { l: "Rapoport β", v: "−0.727" },
                  ].map(m => (
                    <div key={m.l} style={{ borderLeft: `2px solid ${T.teal}20`, paddingLeft: 12 }}>
                      <p style={{ fontFamily: T.mono, color: T.text2, fontSize: 10, letterSpacing: "0.1em", marginBottom: 4 }}>
                        {m.l.toUpperCase()}
                      </p>
                      <p style={{ fontFamily: T.mono, color: T.teal, fontSize: 16, fontWeight: 600 }}>{m.v}</p>
                    </div>
                  ))}
                </div>
                <div className="flex flex-col gap-3">
                  <a href="#" style={{
                    display: "block", textAlign: "center", padding: "12px",
                    background: T.teal, color: T.bg0,
                    fontFamily: T.mono, fontSize: 12, letterSpacing: "0.08em",
                  }}>
                    READ FULL PAPER →
                  </a>
                  <a href="#" style={{
                    display: "block", textAlign: "center", padding: "12px",
                    border: `1px solid ${T.border}`, color: T.text1,
                    fontFamily: T.mono, fontSize: 12, letterSpacing: "0.08em",
                  }}>
                    CITE THIS WORK
                  </a>
                </div>
              </div>
            </div>
          </div>
        </Section>
      </div>

      {/* ── FOOTER ────────────────────────────────────────────────────── */}
      <footer style={{
        background: T.bgDark, borderTop: `1px solid ${T.borderDark}`,
        padding: "24px", textAlign: "center",
      }}>
        <p style={{ fontFamily: T.mono, color: T.dtxt2, fontSize: 11, letterSpacing: "0.08em", lineHeight: 1.8 }}>
          © 2025 · MIGRATORY BUTTERFLY MACROECOLOGY STUDY<br />
          DATA: LEPTRAITS 1.1.0 · WORLDCLIM 2.1 · CLGS · KAWAHARA ET AL. 2023
        </p>
      </footer>
    </div>
  );
}
