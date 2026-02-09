"""Apple-style web frontend for uht-tooling (NiceGUI)."""

from __future__ import annotations

import contextlib
import logging
import socket
from typing import Optional

from nicegui import ui

from uht_tooling.web.app import register_routes

_LOGGER = logging.getLogger("uht_tooling.web")


def _port_is_available(host: str, port: int) -> bool:
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            sock.bind((host, port))
        except OSError:
            return False
    return True


def _find_port(host: str, preferred: Optional[int]) -> int:
    if preferred is not None and _port_is_available(host, preferred):
        return preferred
    if preferred is not None:
        _LOGGER.warning(
            "Port %s unavailable on %s; searching for an open port.", preferred, host
        )
    with contextlib.closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
        sock.bind((host, 0))
        return sock.getsockname()[1]


def launch_web_gui(
    host: str = "127.0.0.1",
    port: Optional[int] = 7860,
) -> None:
    """Start the NiceGUI-based Apple-style GUI."""
    resolved_port = _find_port(host, port)
    _LOGGER.info("Starting uht-tooling web GUI on http://%s:%s", host, resolved_port)
    register_routes()
    ui.run(
        host=host,
        port=resolved_port,
        title="uht-tooling",
        reload=False,
        show=True,
    )
