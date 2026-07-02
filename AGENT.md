# Agent Instructions

This repository contains Quarto documents, R scripts, model outputs, and supporting data for CCSBT southern bluefin tuna assessment work.

## Working Guidelines

- Prefer small, focused changes that preserve the existing project layout.
- Check the relevant `.qmd`, `.R`, and data files before changing analysis behavior.
- Do not delete or regenerate large outputs, cache directories, PDFs, HTML files, `.rda`, `.rds`, or model artifacts unless the user explicitly asks.
- Treat existing uncommitted changes as user work. Do not revert them unless requested.
- Keep generated-document edits separate from source edits when possible; update source `.qmd` or `.R` files first.

## Verification

- For Quarto changes, use `quarto render <file>` when feasible.
- For R script changes, run the narrowest relevant R command or script section first.
- If a full model run or document render is too expensive, state what was and was not verified.

## Style

- Follow the style of nearby R and Quarto code.
- Use clear variable names and avoid broad refactors unless they directly support the requested change.
- Keep comments concise and focused on analysis assumptions, data provenance, or non-obvious modeling logic.
