"""Apple-inspired design system: CSS custom properties, glass effects, animations."""

APPLE_CSS = """
/* ── Typography ────────────────────────────────────────────────── */
:root {
    --font-sans: -apple-system, BlinkMacSystemFont, "SF Pro Display",
                 "Helvetica Neue", Helvetica, Arial, sans-serif;
    --font-mono: "SF Mono", SFMono-Regular, ui-monospace, Menlo, monospace;
}

/* ── Light-mode tokens (default) ───────────────────────────────── */
:root {
    --bg-primary:   #ffffff;
    --bg-secondary: #f5f5f7;
    --bg-glass:     rgba(255, 255, 255, 0.72);
    --bg-glass-heavy: rgba(255, 255, 255, 0.85);
    --text-primary:   #1d1d1f;
    --text-secondary: #6e6e73;
    --accent:       #0071e3;
    --accent-hover: #0077ed;
    --border:       rgba(0, 0, 0, 0.08);
    --border-input: rgba(0, 0, 0, 0.12);
    --shadow-sm:    0 1px 3px rgba(0, 0, 0, 0.04);
    --shadow-md:    0 4px 14px rgba(0, 0, 0, 0.06);
    --shadow-lg:    0 12px 40px rgba(0, 0, 0, 0.12);
    --radius-input: 8px;
    --radius-card:  12px;
    --radius-hero:  20px;
    --radius-pill:  980px;
    --sidebar-width: 260px;
    --content-max-width: 820px;
}

/* ── Dark-mode tokens ──────────────────────────────────────────── */
[data-theme="dark"] {
    --bg-primary:   #1d1d1f;
    --bg-secondary: #2d2d2f;
    --bg-glass:     rgba(29, 29, 31, 0.72);
    --bg-glass-heavy: rgba(29, 29, 31, 0.85);
    --text-primary:   #f5f5f7;
    --text-secondary: #a1a1a6;
    --accent:       #2997ff;
    --accent-hover: #40a9ff;
    --border:       rgba(255, 255, 255, 0.08);
    --border-input: rgba(255, 255, 255, 0.12);
    --shadow-sm:    0 1px 3px rgba(0, 0, 0, 0.2);
    --shadow-md:    0 4px 14px rgba(0, 0, 0, 0.25);
    --shadow-lg:    0 12px 40px rgba(0, 0, 0, 0.4);
}

/* ── Animations ────────────────────────────────────────────────── */
@keyframes slideUpFadeIn {
    from { opacity: 0; transform: translateY(16px); }
    to   { opacity: 1; transform: translateY(0); }
}

@keyframes shimmer {
    0%   { background-position: -200% 0; }
    100% { background-position: 200% 0; }
}

/* ── Base reset & typography ───────────────────────────────────── */
body, .nicegui-content {
    font-family: var(--font-sans) !important;
    background: var(--bg-secondary) !important;
    color: var(--text-primary) !important;
    margin: 0;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
}

/* ── Sidebar ───────────────────────────────────────────────────── */
.apple-sidebar {
    position: fixed;
    top: 0;
    left: 0;
    width: var(--sidebar-width);
    height: 100vh;
    background: var(--bg-glass-heavy);
    backdrop-filter: blur(20px) saturate(180%);
    -webkit-backdrop-filter: blur(20px) saturate(180%);
    border-right: 1px solid var(--border);
    z-index: 100;
    display: flex;
    flex-direction: column;
    overflow-y: auto;
    padding: 24px 16px 16px;
    box-sizing: border-box;
}

.apple-sidebar .sidebar-brand {
    font-size: 20px;
    font-weight: 700;
    letter-spacing: -0.4px;
    color: var(--text-primary);
    margin-bottom: 6px;
}

.apple-sidebar .sidebar-subtitle {
    font-size: 12px;
    color: var(--text-secondary);
    margin-bottom: 28px;
}

.apple-sidebar .sidebar-section {
    font-size: 11px;
    font-weight: 600;
    letter-spacing: 0.5px;
    text-transform: uppercase;
    color: var(--text-secondary);
    margin: 20px 0 8px 8px;
}

.apple-sidebar .sidebar-item {
    display: flex;
    align-items: center;
    padding: 8px 12px;
    border-radius: 8px;
    font-size: 14px;
    font-weight: 500;
    color: var(--text-primary);
    cursor: pointer;
    transition: background 0.15s ease, color 0.15s ease;
    text-decoration: none;
    margin-bottom: 2px;
    user-select: none;
}

.apple-sidebar .sidebar-item .sidebar-icon,
.apple-sidebar .sidebar-item .q-icon,
.apple-sidebar .sidebar-item .material-icons,
.apple-sidebar .sidebar-item i {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    flex: 0 0 18px;
    width: 18px;
    height: 18px;
    line-height: 18px;
    margin-right: 10px;
}

.apple-sidebar .sidebar-item .sidebar-label {
    line-height: 18px;
}

.apple-sidebar .sidebar-item:hover {
    background: var(--border);
}

.apple-sidebar .sidebar-item.active {
    background: var(--accent);
    color: #ffffff;
}

.sidebar-spacer {
    flex: 1;
}

.theme-toggle {
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 8px;
    padding: 8px 12px;
    border-radius: 8px;
    cursor: pointer;
    font-size: 13px;
    color: var(--text-secondary);
    transition: background 0.15s ease;
}

.theme-toggle:hover {
    background: var(--border);
}

/* ── Content area ──────────────────────────────────────────────── */
.apple-content {
    margin-left: var(--sidebar-width);
    padding: 40px 48px;
    min-height: 100vh;
    box-sizing: border-box;
}

.apple-content-inner {
    max-width: var(--content-max-width);
    margin: 0 auto;
}

/* ── Card ──────────────────────────────────────────────────────── */
.apple-card {
    background: var(--bg-glass);
    backdrop-filter: blur(20px) saturate(180%);
    -webkit-backdrop-filter: blur(20px) saturate(180%);
    border: 1px solid var(--border);
    border-radius: var(--radius-card);
    box-shadow: var(--shadow-md);
    padding: 28px 32px;
    margin-bottom: 24px;
    animation: slideUpFadeIn 0.45s ease-out both;
}

.apple-card.delay-1 { animation-delay: 0.08s; }
.apple-card.delay-2 { animation-delay: 0.16s; }
.apple-card.delay-3 { animation-delay: 0.24s; }

.apple-card-title {
    font-size: 22px;
    font-weight: 700;
    letter-spacing: -0.3px;
    color: var(--text-primary);
    margin: 0 0 4px 0;
}

.apple-card-subtitle {
    font-size: 14px;
    color: var(--text-secondary);
    margin: 0 0 20px 0;
    line-height: 1.5;
}

/* ── Form inputs ───────────────────────────────────────────────── */
.apple-field label,
.apple-field .q-field__label {
    font-size: 13px !important;
    font-weight: 600 !important;
    color: var(--text-secondary) !important;
    letter-spacing: 0.1px;
}

.apple-field .q-field__control {
    border: 1px solid var(--border-input) !important;
    border-radius: var(--radius-input) !important;
    background: var(--bg-primary) !important;
    transition: border-color 0.2s ease, box-shadow 0.2s ease;
    font-size: 15px !important;
}

.apple-field .q-field__control:focus-within {
    border-color: var(--accent) !important;
    box-shadow: 0 0 0 3px rgba(0, 113, 227, 0.15) !important;
}

[data-theme="dark"] .apple-field .q-field__control:focus-within {
    box-shadow: 0 0 0 3px rgba(41, 151, 255, 0.2) !important;
}

.apple-field .q-field__control::before,
.apple-field .q-field__control::after {
    border: none !important;
}

.apple-field .q-field__native,
.apple-field textarea {
    color: var(--text-primary) !important;
    font-family: var(--font-sans) !important;
    font-size: 15px !important;
}

.apple-field .q-field__native::placeholder,
.apple-field textarea::placeholder {
    color: var(--text-secondary) !important;
    opacity: 0.5;
}

/* Hide placeholders when label is inline to prevent overlap */
.apple-field:not(.q-field--focused) .q-field__native::placeholder,
.apple-field:not(.q-field--focused) textarea::placeholder {
    opacity: 0;
}

/* ── Number field arrows ───────────────────────────────────────── */
.apple-field .q-field__append {
    color: var(--text-secondary) !important;
}

/* ── Slider ────────────────────────────────────────────────────── */
.apple-slider .q-slider__track-container {
    height: 4px !important;
}

.apple-slider .q-slider__track {
    border-radius: 2px;
}

.apple-slider .q-slider__thumb {
    width: 20px !important;
    height: 20px !important;
    background: var(--bg-primary) !important;
    border: 2px solid var(--accent) !important;
    box-shadow: var(--shadow-md) !important;
}

.apple-slider .q-slider__thumb::after {
    display: none;
}

.apple-slider .q-slider__thumb-shape {
    color: var(--accent) !important;
}

/* ── Buttons ───────────────────────────────────────────────────── */
.apple-btn {
    background: var(--accent) !important;
    color: #ffffff !important;
    border: none !important;
    border-radius: var(--radius-pill) !important;
    padding: 10px 28px !important;
    font-size: 15px !important;
    font-weight: 600 !important;
    letter-spacing: -0.1px !important;
    text-transform: none !important;
    box-shadow: var(--shadow-sm) !important;
    transition: transform 0.12s ease, box-shadow 0.12s ease, background 0.15s ease !important;
    cursor: pointer;
}

.apple-btn:hover {
    background: var(--accent-hover) !important;
    box-shadow: var(--shadow-md) !important;
}

.apple-btn:active {
    transform: scale(0.98) !important;
}

.apple-btn-secondary {
    background: var(--bg-primary) !important;
    color: var(--accent) !important;
    border: 1px solid var(--border-input) !important;
    border-radius: var(--radius-pill) !important;
    padding: 8px 22px !important;
    font-size: 14px !important;
    font-weight: 600 !important;
    text-transform: none !important;
    box-shadow: none !important;
    transition: transform 0.12s ease, background 0.15s ease !important;
    cursor: pointer;
}

.apple-btn-secondary:hover {
    background: var(--bg-secondary) !important;
}

.apple-btn-secondary:active {
    transform: scale(0.98) !important;
}

/* ── Upload dropzone ───────────────────────────────────────────── */
.apple-upload .q-uploader {
    border: 2px dashed var(--border-input) !important;
    border-radius: var(--radius-card) !important;
    background: var(--bg-primary) !important;
    transition: border-color 0.2s ease, background 0.2s ease;
    width: 100%;
    box-shadow: none !important;
}

.apple-upload .q-uploader:hover {
    border-color: var(--accent) !important;
    background: rgba(0, 113, 227, 0.03) !important;
}

[data-theme="dark"] .apple-upload .q-uploader:hover {
    background: rgba(41, 151, 255, 0.05) !important;
}

.apple-upload .q-uploader__header {
    background: transparent !important;
    color: var(--text-secondary) !important;
}

.apple-upload .q-uploader__title,
.apple-upload .q-uploader__subtitle {
    color: var(--text-secondary) !important;
}

/* ── Table ─────────────────────────────────────────────────────── */
.apple-table .q-table {
    background: transparent !important;
    border: 1px solid var(--border) !important;
    border-radius: var(--radius-input) !important;
    overflow: hidden;
}

.apple-table .q-table thead tr th {
    font-size: 12px !important;
    font-weight: 600 !important;
    text-transform: uppercase !important;
    letter-spacing: 0.5px !important;
    color: var(--text-secondary) !important;
    background: var(--bg-secondary) !important;
    border-bottom: 1px solid var(--border) !important;
}

.apple-table .q-table tbody tr td {
    font-size: 14px !important;
    color: var(--text-primary) !important;
    border-bottom: 1px solid var(--border) !important;
}

.apple-table .q-table tbody tr:hover td {
    background: rgba(0, 113, 227, 0.03) !important;
}

[data-theme="dark"] .apple-table .q-table tbody tr:hover td {
    background: rgba(41, 151, 255, 0.05) !important;
}

/* ── Markdown ──────────────────────────────────────────────────── */
.apple-markdown {
    font-family: var(--font-sans);
    font-size: 15px;
    line-height: 1.65;
    color: var(--text-primary);
}

.apple-markdown h1, .apple-markdown h2, .apple-markdown h3 {
    font-weight: 700;
    letter-spacing: -0.3px;
    color: var(--text-primary);
}

.apple-markdown code {
    font-family: var(--font-mono);
    font-size: 13px;
    background: var(--bg-secondary);
    padding: 2px 6px;
    border-radius: 4px;
}

.apple-markdown pre {
    background: var(--bg-secondary);
    border: 1px solid var(--border);
    border-radius: var(--radius-input);
    padding: 16px;
    overflow-x: auto;
}

.apple-markdown pre code {
    background: none;
    padding: 0;
}

.apple-markdown table {
    border-collapse: collapse;
    width: 100%;
    font-size: 14px;
}

.apple-markdown table th,
.apple-markdown table td {
    border: 1px solid var(--border);
    padding: 8px 12px;
    text-align: left;
}

.apple-markdown table th {
    background: var(--bg-secondary);
    font-weight: 600;
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 0.4px;
    color: var(--text-secondary);
}

/* ── Progress shimmer bar ──────────────────────────────────────── */
.apple-progress-bar {
    height: 3px;
    width: 100%;
    border-radius: 2px;
    background: linear-gradient(
        90deg,
        var(--bg-secondary) 0%,
        var(--accent) 50%,
        var(--bg-secondary) 100%
    );
    background-size: 200% 100%;
    animation: shimmer 1.8s infinite linear;
    margin: 12px 0;
}

/* ── Error banner ──────────────────────────────────────────────── */
.apple-error {
    background: rgba(255, 59, 48, 0.08);
    border: 1px solid rgba(255, 59, 48, 0.2);
    border-radius: var(--radius-input);
    padding: 12px 16px;
    color: #ff3b30;
    font-size: 14px;
    font-weight: 500;
    margin: 12px 0;
}

[data-theme="dark"] .apple-error {
    background: rgba(255, 69, 58, 0.12);
    border-color: rgba(255, 69, 58, 0.25);
    color: #ff453a;
}

/* ── Quasar overrides for Apple feel ───────────────────────────── */
.q-field--outlined .q-field__control::before {
    border: none !important;
}
.q-field--outlined .q-field__control::after {
    border: none !important;
}

/* Remove Quasar ripple effect for cleaner look */
.q-ripple {
    display: none !important;
}

/* ── Responsive ────────────────────────────────────────────────── */
@media (max-width: 860px) {
    .apple-sidebar {
        width: 100%;
        height: auto;
        position: relative;
        border-right: none;
        border-bottom: 1px solid var(--border);
    }
    .apple-content {
        margin-left: 0;
        padding: 24px 20px;
    }
}
"""
