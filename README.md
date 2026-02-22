# Spring Design Dashboard

A Streamlit dashboard for analyzing and optimizing compression spring designs, with a focus on blaster applications.

## Features

### Analysis Page
- Input spring parameters (OD, free length, wire diameter, active coils, end type, wire type)
- Spring presets for common springs (SFX3–SFX8, 5kg LS, 8kg LS)
- Calculated specs: spring rate, inner diameter, total coils, spring index, pitch, weight, and more
- Safe/solid load table with WB Jones 45% safe rule assessment
- Detailed calculations with rendered LaTeX formulas
- Force–deflection diagram
- Operating-point stress check
- FPS estimation and efficiency back-calculation with blaster presets (BABP, Hijinx, Lynx, Minx, Lonx)
- All parameters encoded in the URL for bookmarking and sharing

### Optimization Page
- Target a specific FPS (at a given efficiency and dart weight) or a spring rate
- Blaster geometry constraints with presets (compress from/to)
- OD constraint (fixed or range)
- Free length (exact or derived from compress-from with settling allowance)
- Brute-force search across all wire diameters and OD values
- Results ranked by safe-to-solid preference and solid height margin
- Pareto-optimal filtering to remove dominated designs
- One-click "Open in Analysis" to inspect any candidate in detail

## Quick Start

```bash
pip install -r requirements.txt
streamlit run spring_dashboard.py
```

## Project Structure

```
spring_dashboard.py    # Entry point — sets up multi-page navigation
spring_helpers.py      # Shared constants, physics functions, UI helpers
pages/
  analysis.py          # Spring analysis page
  optimization.py      # Spring optimization/search page
```
