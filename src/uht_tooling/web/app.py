"""App shell: sidebar navigation, routing, dark-mode toggle."""

from __future__ import annotations

from nicegui import app, ui

from uht_tooling.web.theme import APPLE_CSS
from uht_tooling.web.pages import (
    nextera,
    slim,
    kld,
    gibson,
    mutation_caller,
    umi_hunter,
    profile_inserts,
    ep_library,
)

# Map of route suffix → (label, page builder)
_NAV_GROUPS: list[tuple[str, list[tuple[str, str, str]]]] = [
    (
        "Primer Design",
        [
            ("/", "Nextera XT", "science"),
            ("/slim", "SLIM", "biotech"),
            ("/kld", "KLD", "loop"),
            ("/gibson", "Gibson", "merge_type"),
        ],
    ),
    (
        "Sequencing Analysis",
        [
            ("/mutation-caller", "Mutation Caller", "troubleshoot"),
            ("/umi-hunter", "UMI Hunter", "qr_code"),
            ("/profile-inserts", "Profile Inserts", "insert_chart"),
            ("/ep-library", "EP Library Profile", "auto_graph"),
        ],
    ),
]


def _build_sidebar(current_path: str) -> None:
    """Render the fixed left sidebar with navigation links."""
    with ui.element("nav").classes("apple-sidebar"):
        ui.label("uht-tooling").classes("sidebar-brand")
        ui.label("Primer design & sequencing analysis").classes("sidebar-subtitle")

        for group_label, items in _NAV_GROUPS:
            ui.label(group_label).classes("sidebar-section")
            for route, label, icon_name in items:
                active = " active" if current_path == route else ""
                with ui.link(target=route).classes(
                    f"sidebar-item{active}"
                ).style("text-decoration: none;"):
                    ui.icon(icon_name).classes("sidebar-icon").style("font-size: 18px;")
                    ui.label(label).classes("sidebar-label")

        ui.element("div").classes("sidebar-spacer")

        # Dark-mode toggle
        dark = ui.dark_mode()

        def _get_saved_theme() -> str:
            try:
                return app.storage.user.get("theme", "light")
            except RuntimeError:
                return "light"

        def _set_saved_theme(theme: str) -> None:
            try:
                app.storage.user["theme"] = theme
            except RuntimeError:
                return

        def apply_saved_theme() -> None:
            theme = _get_saved_theme()
            if theme == "dark":
                dark.enable()
            else:
                dark.disable()
            ui.run_javascript(
                f"document.documentElement.setAttribute('data-theme', '{theme}')"
            )

        apply_saved_theme()

        def toggle_theme():
            current = _get_saved_theme()
            new_theme = "dark" if current != "dark" else "light"
            _set_saved_theme(new_theme)
            if new_theme == "dark":
                dark.enable()
            else:
                dark.disable()
            ui.run_javascript(
                f"document.documentElement.setAttribute('data-theme', '{new_theme}')"
            )

        with ui.element("div").classes("theme-toggle").on("click", toggle_theme):
            ui.icon("dark_mode").style("font-size: 18px;")
            ui.label("Toggle Dark Mode")


def _page_wrapper(current_path: str):
    """Decorator-style context: injects CSS, sidebar, content container."""
    ui.add_head_html(f"<style>{APPLE_CSS}</style>")
    _build_sidebar(current_path)
    content = ui.element("main").classes("apple-content")
    return content


# ── Route registration ─────────────────────────────────────────────────────


def register_routes() -> None:
    """Register all NiceGUI page routes."""

    @ui.page("/")
    async def page_nextera():
        with _page_wrapper("/"):
            with ui.element("div").classes("apple-content-inner"):
                await nextera.render()

    @ui.page("/slim")
    async def page_slim():
        with _page_wrapper("/slim"):
            with ui.element("div").classes("apple-content-inner"):
                await slim.render()

    @ui.page("/kld")
    async def page_kld():
        with _page_wrapper("/kld"):
            with ui.element("div").classes("apple-content-inner"):
                await kld.render()

    @ui.page("/gibson")
    async def page_gibson():
        with _page_wrapper("/gibson"):
            with ui.element("div").classes("apple-content-inner"):
                await gibson.render()

    @ui.page("/mutation-caller")
    async def page_mutation_caller():
        with _page_wrapper("/mutation-caller"):
            with ui.element("div").classes("apple-content-inner"):
                await mutation_caller.render()

    @ui.page("/umi-hunter")
    async def page_umi_hunter():
        with _page_wrapper("/umi-hunter"):
            with ui.element("div").classes("apple-content-inner"):
                await umi_hunter.render()

    @ui.page("/profile-inserts")
    async def page_profile_inserts():
        with _page_wrapper("/profile-inserts"):
            with ui.element("div").classes("apple-content-inner"):
                await profile_inserts.render()

    @ui.page("/ep-library")
    async def page_ep_library():
        with _page_wrapper("/ep-library"):
            with ui.element("div").classes("apple-content-inner"):
                await ep_library.render()
