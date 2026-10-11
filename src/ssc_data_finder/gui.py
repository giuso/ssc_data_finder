"""Entry point for `ssc-data-finder-gui`: launches the Streamlit app and opens it in the browser."""
from __future__ import annotations

import sys
import threading
import webbrowser
from pathlib import Path

DEFAULT_PORT = 8501


def main() -> None:
    from streamlit.web import cli as stcli

    args = sys.argv[1:]
    port = next((a.split("=", 1)[1] for a in args if a.startswith("--server.port=")), DEFAULT_PORT)
    # Headless skips Streamlit's first-run email prompt; open the browser ourselves instead.
    threading.Timer(2.0, webbrowser.open, args=[f"http://localhost:{port}"]).start()
    sys.argv = [
        "streamlit",
        "run",
        str(Path(__file__).with_name("app.py")),
        "--server.headless=true",
        "--browser.gatherUsageStats=false",
        *args,
    ]
    sys.exit(stcli.main())


if __name__ == "__main__":
    main()
