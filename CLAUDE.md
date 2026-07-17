# CLAUDE.md

Guidance for Claude when working in this research repo.

## Topics

- Write each topic as a markdown file under `topics/<topic-name>/README.md`.
- Keep it **short**: bullet points, quick read. Only go long-form if I explicitly ask for depth.
- Code samples are **Python**.
- Add graphs/plots when I ask for a visualisation (matplotlib is the default).
- For video / animated examples, use [Manim](https://www.manim.community/).
- Whenever you add, remove, or rename a topic, **update the root `README.md`** so its "Topics" section lists every topic in the `topics/` directory (link + one-line description).

## Project

- This is a **uv** project. Manage dependencies with `uv add` / `uv add --dev`, run code with `uv run`.
- **ruff** runs via pre-commit. Keep code clean so hooks pass.

## Git

- Committing straight to `main` is fine.
- **Never erase history** — no force-pushes, no rewriting past commits.
