from __future__ import annotations

import sys
import webbrowser
from pathlib import Path


def _find_index_html() -> Path:
    """Return the local index.html path for an editable repository install."""
    candidates = [
        Path.cwd() / "index.html",
        Path(__file__).resolve().parent.parent / "index.html",
    ]

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        "index.html was not found. Run this command from the 4MRNA-G4 repository "
        "or install the repository with `python -m pip install -e .`."
    )


def main() -> None:
    try:
        index_html = _find_index_html()
    except FileNotFoundError as exc:
        print(f"4MRNA-G4: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    webbrowser.open(index_html.resolve().as_uri())


if __name__ == "__main__":
    main()
