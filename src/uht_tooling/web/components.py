"""Reusable Apple-styled component wrappers around NiceGUI elements."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Callable, List, Optional

from nicegui import ui


# ── Cards ──────────────────────────────────────────────────────────────────


@contextmanager
def apple_card(title: str, subtitle: str = "", delay: int = 0):
    """Frosted-glass card with optional staggered animation delay."""
    delay_cls = f" delay-{delay}" if delay else ""
    with ui.element("div").classes(f"apple-card{delay_cls}") as card:
        if title:
            ui.label(title).classes("apple-card-title")
        if subtitle:
            ui.label(subtitle).classes("apple-card-subtitle")
        yield card


# ── Text inputs ────────────────────────────────────────────────────────────


def apple_input(
    label: str,
    placeholder: str = "",
    *,
    value: str = "",
) -> ui.input:
    inp = ui.input(label=label, placeholder=placeholder, value=value)
    inp.classes("apple-field w-full")
    inp.props('outlined dense')
    return inp


def apple_textarea(
    label: str,
    placeholder: str = "",
    *,
    value: str = "",
    rows: int = 4,
) -> ui.textarea:
    ta = ui.textarea(label=label, placeholder=placeholder, value=value)
    ta.classes("apple-field w-full")
    ta.props(f'outlined dense rows={rows}')
    return ta


# ── Number fields ──────────────────────────────────────────────────────────


def apple_number(
    label: str,
    value: float = 0,
    *,
    min_val: Optional[float] = None,
    max_val: Optional[float] = None,
    step: float = 1,
    precision: int = 0,
) -> ui.number:
    kwargs: dict[str, Any] = {"label": label, "value": value, "step": step}
    if min_val is not None:
        kwargs["min"] = min_val
    if max_val is not None:
        kwargs["max"] = max_val
    if precision == 0:
        kwargs["format"] = "%.0f"
    n = ui.number(**kwargs)
    n.classes("apple-field w-full")
    n.props('outlined dense')
    return n


# ── Sliders ────────────────────────────────────────────────────────────────


def apple_slider(
    label: str,
    min_val: float,
    max_val: float,
    value: float,
    step: float = 0.01,
) -> ui.slider:
    ui.label(label).classes("text-xs font-semibold").style(
        "color: var(--text-secondary); margin-bottom: 2px;"
    )
    s = ui.slider(min=min_val, max=max_val, value=value, step=step)
    s.classes("apple-slider w-full")
    return s


# ── Upload ─────────────────────────────────────────────────────────────────


def apple_upload(
    label: str,
    *,
    extensions: Optional[List[str]] = None,
    multiple: bool = False,
    on_upload: Optional[Callable] = None,
    auto_upload: bool = True,
) -> ui.upload:
    """Dashed-border dropzone upload."""
    ui.label(label).classes("text-xs font-semibold").style(
        "color: var(--text-secondary); margin-bottom: 4px; margin-top: 8px;"
    )
    u = ui.upload(
        multiple=multiple,
        auto_upload=auto_upload,
        on_upload=on_upload,
    )
    u.classes("apple-upload w-full")
    return u


# ── Table ──────────────────────────────────────────────────────────────────


def apple_table(
    columns: list[dict],
    rows: list[dict],
    *,
    row_key: str = "name",
) -> ui.table:
    t = ui.table(columns=columns, rows=rows, row_key=row_key)
    t.classes("apple-table w-full")
    return t


# ── Buttons ────────────────────────────────────────────────────────────────


def apple_button(
    label: str,
    *,
    on_click: Optional[Callable] = None,
) -> ui.button:
    b = ui.button(label, on_click=on_click)
    b.classes("apple-btn")
    b.props('unelevated no-caps')
    return b


def apple_button_secondary(
    label: str,
    *,
    on_click: Optional[Callable] = None,
) -> ui.button:
    b = ui.button(label, on_click=on_click)
    b.classes("apple-btn-secondary")
    b.props('unelevated no-caps')
    return b


# ── Markdown ───────────────────────────────────────────────────────────────


def apple_markdown(content: str = "") -> ui.markdown:
    m = ui.markdown(content)
    m.classes("apple-markdown")
    return m


# ── Progress shimmer ──────────────────────────────────────────────────────


def apple_progress() -> ui.element:
    el = ui.element("div").classes("apple-progress-bar")
    el.set_visibility(False)
    return el


# ── Download helper ────────────────────────────────────────────────────────


def apple_download_button(
    label: str = "Download Results",
    zip_path: str = "",
) -> ui.button:
    btn = ui.button(label, on_click=lambda: ui.download(zip_path))
    btn.classes("apple-btn-secondary")
    btn.props('unelevated no-caps')
    btn.set_visibility(False)
    return btn
