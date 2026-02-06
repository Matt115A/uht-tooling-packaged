"""Global configuration file support."""
import os
from pathlib import Path
from typing import Any, Dict, Optional

try:
    import yaml

    HAVE_YAML = True
except ImportError:
    HAVE_YAML = False


DEFAULT_CONFIG_PATHS = [
    Path.home() / ".uht-tooling.yaml",
    Path.home() / ".config" / "uht-tooling" / "config.yaml",
    Path(".uht-tooling.yaml"),
]


def find_config_file() -> Optional[Path]:
    """
    Find a configuration file from environment variable or default locations.

    Search order:
    1. $UHT_TOOLING_CONFIG environment variable
    2. ~/.uht-tooling.yaml
    3. ~/.config/uht-tooling/config.yaml
    4. .uht-tooling.yaml (current directory)

    Returns:
        Path to the config file if found, None otherwise.
    """
    # Check environment variable first
    env_path = os.environ.get("UHT_TOOLING_CONFIG")
    if env_path:
        path = Path(env_path)
        if path.exists():
            return path

    # Check default locations
    for path in DEFAULT_CONFIG_PATHS:
        if path.exists():
            return path

    return None


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load YAML configuration, auto-discovering if path not provided.

    Args:
        config_path: Explicit path to config file. If None, auto-discover.

    Returns:
        Dictionary containing configuration. Empty dict if no config found
        or if YAML is not available.
    """
    if not HAVE_YAML:
        return {}

    if config_path is None:
        config_path = find_config_file()

    if config_path is not None:
        config_path = Path(config_path)

    if config_path is None or not config_path.exists():
        return {}

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
            return config if isinstance(config, dict) else {}
    except Exception:
        return {}


def get_option(
    config: Dict[str, Any],
    key: str,
    cli_value: Any,
    default: Any = None,
    workflow: Optional[str] = None,
) -> Any:
    """
    Get an option with precedence: CLI > workflow-specific config > global config > default.

    Args:
        config: Configuration dictionary from load_config().
        key: The option key to look up.
        cli_value: Value from CLI (takes precedence if not None).
        default: Default value if not found anywhere.
        workflow: Optional workflow name for workflow-specific defaults.

    Returns:
        The resolved option value.
    """
    # CLI value always takes precedence if explicitly provided
    if cli_value is not None:
        return cli_value

    # Check workflow-specific defaults
    if workflow:
        workflow_defaults = config.get("defaults", {}).get(workflow, {})
        if key in workflow_defaults:
            return workflow_defaults[key]

    # Check global paths config
    paths_config = config.get("paths", {})
    if key in paths_config:
        value = paths_config[key]
        # Expand ~ in paths
        if isinstance(value, str):
            return os.path.expanduser(value)
        return value

    # Check top-level config
    if key in config:
        return config[key]

    return default


def get_workflow_defaults(config: Dict[str, Any], workflow: str) -> Dict[str, Any]:
    """
    Get all default values for a specific workflow.

    Args:
        config: Configuration dictionary from load_config().
        workflow: Workflow name (e.g., "mutation_caller", "umi_hunter").

    Returns:
        Dictionary of default values for the workflow.
    """
    return config.get("defaults", {}).get(workflow, {})
