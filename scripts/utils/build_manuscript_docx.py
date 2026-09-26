#!/usr/bin/env python3
"""
Convert manuscript/manuscript.md to manuscript/manuscript.docx.

Deliberate markdown subset parser (not a general-purpose converter):
headings (#, ##, ###), blockquote notes, bold (**)/italic (*), bullet and
numbered lists, horizontal rules, pipe tables (rendered as real Word tables),
and image references (![alt](path), embedded as real pictures). Good enough
for this one document.

Usage:
    python3 build_manuscript_docx.py [src.md] [dst.docx]
Defaults to manuscript/manuscript.md -> manuscript/manuscript.docx (paths
resolved relative to this script's location, so it runs from any cwd).
Image paths in the markdown are resolved relative to the markdown file's
own directory.
"""
import os
import re
import sys

from docx import Document
from docx.shared import Pt, RGBColor, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_SRC = os.path.join(REPO_ROOT, "manuscript", "manuscript.md")

SRC = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SRC
DST = sys.argv[2] if len(sys.argv) > 2 else SRC.replace(".md", ".docx")
SRC_DIR = os.path.dirname(os.path.abspath(SRC))

IMG_RE = re.compile(r'^!\[([^\]]*)\]\(([^)]+)\)$')
MAX_IMG_WIDTH_IN = 6.3


def add_runs_with_bold(paragraph, text):
    """Split on **bold** and *italic* markers and add runs accordingly."""
    parts = re.split(r'(\*\*.*?\*\*|\*[^*]+?\*)', text)
    for part in parts:
        if not part:
            continue
        if part.startswith('**') and part.endswith('**'):
            run = paragraph.add_run(part[2:-2])
            run.bold = True
        elif part.startswith('*') and part.endswith('*') and len(part) > 1:
            run = paragraph.add_run(part[1:-1])
            run.italic = True
        else:
            paragraph.add_run(part)


def set_cell_shading(cell, hex_color):
    shd = OxmlElement('w:shd')
    shd.set(qn('w:fill'), hex_color)
    cell._tc.get_or_add_tcPr().append(shd)


def parse_table_block(lines, start_idx):
    """Given lines[start_idx] is the header row of a pipe table, consume the
    header, the '---' separator, and all following data rows. Returns
    (rows, next_idx) where rows is a list of lists of cell-text strings."""
    def split_row(line):
        cells = line.strip().strip('|').split('|')
        return [c.strip() for c in cells]

    rows = [split_row(lines[start_idx])]
    idx = start_idx + 1
    if idx < len(lines) and set(lines[idx].strip()) <= set('|:- '):
        idx += 1
    while idx < len(lines) and lines[idx].strip().startswith('|'):
        rows.append(split_row(lines[idx]))
        idx += 1
    return rows, idx


def add_table(doc, rows):
    n_cols = len(rows[0])
    table = doc.add_table(rows=len(rows), cols=n_cols)
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for r, row_cells in enumerate(rows):
        for c, cell_text in enumerate(row_cells):
            if c >= n_cols:
                continue
            cell = table.cell(r, c)
            cell.paragraphs[0].text = ''
            p = cell.paragraphs[0]
            add_runs_with_bold(p, cell_text)
            if r == 0:
                for run in p.runs:
                    run.bold = True
                set_cell_shading(cell, 'D9E2F3')
            for run in p.runs:
                run.font.size = Pt(9.5)
    doc.add_paragraph()


def add_image(doc, alt_text, rel_path):
    img_path = os.path.normpath(os.path.join(SRC_DIR, rel_path))
    if not os.path.isfile(img_path):
        p = doc.add_paragraph()
        run = p.add_run(f"[MISSING IMAGE: {rel_path}]")
        run.bold = True
        run.font.color.rgb = RGBColor(0xCC, 0x00, 0x00)
        print(f"  WARNING: image not found: {img_path}", file=sys.stderr)
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(img_path, width=Inches(MAX_IMG_WIDTH_IN))
    if alt_text:
        cap = doc.add_paragraph()
        cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
        cap_run = cap.add_run(alt_text)
        cap_run.italic = True
        cap_run.font.size = Pt(9)
        cap_run.font.color.rgb = RGBColor(0x44, 0x44, 0x44)


def main():
    with open(SRC, encoding='utf-8') as f:
        lines = f.read().split('\n')

    doc = Document()
    style = doc.styles['Normal']
    style.font.name = 'Times New Roman'
    style.font.size = Pt(11)

    n_images = 0
    n_tables = 0

    i = 0
    n = len(lines)
    while i < n:
        raw = lines[i]
        line = raw.rstrip()
        if not line.strip():
            i += 1
            continue

        img_match = IMG_RE.match(line.strip())

        if line.startswith('# '):
            h = doc.add_heading(line[2:].strip(), level=0)
            h.alignment = WD_ALIGN_PARAGRAPH.CENTER
        elif line.startswith('## '):
            doc.add_heading(line[3:].strip(), level=1)
        elif line.startswith('### '):
            doc.add_heading(line[4:].strip(), level=2)
        elif line.startswith('>'):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Pt(24)
            run_text = line.lstrip('>').strip()
            add_runs_with_bold(p, run_text)
            for run in p.runs:
                run.italic = True
                run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
        elif line.startswith('---'):
            doc.add_paragraph('_' * 60)
        elif img_match:
            add_image(doc, img_match.group(1), img_match.group(2))
            n_images += 1
        elif line.strip().startswith('|'):
            rows, next_i = parse_table_block(lines, i)
            add_table(doc, rows)
            n_tables += 1
            i = next_i
            continue
        elif re.match(r'^\d+\.\s', line) or line.startswith('- '):
            p = doc.add_paragraph(style='List Bullet' if line.startswith('- ') else 'List Number')
            text = line.split('. ', 1)[-1] if line[0].isdigit() else line[2:]
            add_runs_with_bold(p, text)
        elif line.startswith('*(') and line.endswith(')*'):
            p = doc.add_paragraph()
            run = p.add_run(line.strip('*'))
            run.italic = True
        else:
            p = doc.add_paragraph()
            add_runs_with_bold(p, line)
        i += 1

    doc.save(DST)
    print(f"Saved: {DST}")
    print(f"Embedded {n_images} images and {n_tables} tables.")


if __name__ == '__main__':
    main()
