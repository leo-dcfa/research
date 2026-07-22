# CLAUDE.md

Guidance for Claude when working in this research repo.

## Topics

- Write each topic as a markdown file under `topics/<topic-name>/README.md`.
- Keep it **short**: bullet points, quick read. Only go long-form if I explicitly ask for depth.
- **Equations: always use LaTeX** — `$...$` inline, `$$...$$` (or ```` ```math ```` fences) for display. Never ASCII/unicode math. Inside tables use `\lvert \rvert` etc. instead of raw `|` so the table isn't broken.
- Code samples are **Python**.
- Add graphs/plots when I ask for a visualisation (matplotlib is the default).
- For video / animated examples, use [Manim](https://www.manim.community/).

## Project

- This is a **uv** project. Manage dependencies with `uv add` / `uv add --dev`, run code with `uv run`.
- **ruff** runs via pre-commit. Keep code clean so hooks pass.

## Git

- Committing straight to `main` is fine.
- **Never erase history** — no force-pushes, no rewriting past commits.
