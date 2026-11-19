"""Compile the combined report and retain links from its embedded foundation pages."""
from pathlib import Path
import subprocess

import pymupdf

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "latex"
BUILD = SOURCE / "build"


def main():
    BUILD.mkdir(exist_ok=True)
    for _ in range(2):
        subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
             "-output-directory=build", "main.tex"],
            cwd=SOURCE, check=True,
        )

    # pdfpages preserves the printed material but discards source annotations.
    with pymupdf.open(BUILD / "main.pdf") as report, pymupdf.open(
        SOURCE / "assets/foundation.pdf"
    ) as foundation:
        offset = next(page - 1 for _, title, page in report.get_toc() if title == "Objective")
        for index, page in enumerate(foundation):
            for link in page.get_links():
                if link["kind"] == pymupdf.LINK_URI:
                    target = {"kind": pymupdf.LINK_URI, "uri": link["uri"]}
                elif link.get("page", -1) >= 0:
                    target = {"kind": pymupdf.LINK_GOTO, "page": offset + link["page"],
                              "to": link["to"], "zoom": 0}
                else:
                    continue
                report[offset + index].insert_link({"from": link["from"], **target})
        output = BUILD / "Final_Report.pdf"
        report.save(output, garbage=4, deflate=True)
    output.replace(ROOT / "Final_Report.pdf")


if __name__ == "__main__":
    main()
