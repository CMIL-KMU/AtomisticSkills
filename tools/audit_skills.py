"""Build a complete source inventory without importing scientific engines or MCP servers."""

import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import yaml

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "docs/integrations"


def source_record(path):
    text = path.read_text()
    tree = ast.parse(text)
    functions = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.append(
                dict(
                    name=node.name,
                    line=node.lineno,
                    signature=ast.unparse(node.args),
                    doc=ast.get_docstring(node),
                    calls=sorted(
                        {
                            ast.unparse(n.func)
                            for n in ast.walk(node)
                            if isinstance(n, ast.Call)
                        }
                    ),
                )
            )
    return dict(
        source=str(path.relative_to(ROOT)),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        imports=sorted(
            {
                ast.unparse(n)
                for n in ast.walk(tree)
                if isinstance(n, (ast.Import, ast.ImportFrom))
            }
        ),
        functions=functions,
        cli=[
            ast.unparse(n)
            for n in ast.walk(tree)
            if isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr == "add_argument"
        ],
    )


def main():
    with (DEST / "review.tsv").open() as f:
        review = {r["skill"]: r for r in csv.DictReader(f, delimiter="\t")}
    skills = sorted((ROOT / ".agents/skills").glob("*/SKILL.md")) + sorted(
        (ROOT / "skills").glob("*/SKILL.md")
    )
    assert (
        len(skills) == 130
        and len(review) == 130
        and {p.parent.name for p in skills} == set(review)
    )
    servers = {
        p.stem.removesuffix("_server"): source_record(p)
        for p in (ROOT / "src/mcp_server").glob("*_server.py")
    }
    rows = []
    for path in skills:
        text = path.read_text()
        name = path.parent.name
        r = review[name]
        front = yaml.safe_load(text.split("---", 2)[1])
        refs = sorted(set(re.findall(r"\bmcp_([a-z0-9]+)_([a-zA-Z0-9_]+)", text)))
        entrypoints = []
        for server, method in refs:
            source = servers.get(server)
            matches = (
                [f for f in source["functions"] if f["name"] == method]
                if source
                else []
            )
            entrypoints.append(
                dict(
                    reference="mcp_" + server + "_" + method,
                    source=source["source"] if source else None,
                    source_sha256=source["sha256"] if source else None,
                    availability="defined" if matches else "unresolved-reference",
                    contracts=matches,
                )
            )
        scripts = [
            source_record(p) for p in sorted((path.parent / "scripts").rglob("*.py"))
        ]
        sections = []
        for heading, body in re.findall(
            r"(?m)^(#{1,4} [^\n]+)\n(.*?)(?=^#{1,4} |\Z)", text, re.S
        ):
            if re.search(
                r"input|output|parameter|argument|constraint|requirement|prerequisite|environment",
                heading,
                re.I,
            ):
                sections.append(dict(heading=heading, text=body.strip()))
        row = dict(
            skill=name,
            source=str(path.relative_to(ROOT)),
            source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            description=front.get("description"),
            classification=r["classification"],
            recipe_integration=dict(
                status=r["classification"],
                adapters=r["adapters"].split(",") if r["adapters"] else [],
                remaining_blockers=r["remaining_blockers"],
            ),
            inputs_units_outputs=r["inputs_outputs"],
            provider_contract_sections=sections,
            environments=sorted(set(re.findall(r"\b[a-z0-9]+-agent\b", text))),
            scripts=scripts,
            mcp_entrypoints=entrypoints,
            callable_availability="source-inspected; see each entrypoint",
            installed_dependencies="See environment-availability receipt; metadata is not an import or execution test",
            execution_verification="synthetic adapter tests"
            if r["classification"] == "partial"
            else "pre-existing receipts"
            if r["classification"] == "existing"
            else "not executed in this audit",
            instructions=text,
        )
        rows.append(row)
    result = dict(
        schema="anvil.skill-inventory/v1",
        skill_count=len(rows),
        source_counts={".agents/skills": 129, "skills": 1},
        review_scope="Full instruction text, local script AST signatures/imports/calls and referenced MCP definitions; no blanket runtime success claim",
        skills=rows,
    )
    (DEST / "skill-inventory.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    )
    table = [
        "# AtomisticSkills integration matrix",
        "",
        "All 130 skills are listed. Partial means only the named bounded operation is integrated; it never means the whole tutorial ran. Candidate means remaining implementation work, not an unavailable dependency or a completed integration. Acquisition/literature/review guidance is not represented as fake calculation success.",
        "",
        "The [machine-readable inventory](skill-inventory.json) retains the full instructions, source SHA256, exact CLI options, function signatures/imports/calls, declared environments and unresolved MCP references for every row. Installed package metadata and actual test receipts are separate evidence.",
        "",
        "| Skill and source | Classification / adapters | Inputs, units → outputs | Exact remaining scope or blocker |",
        "|---|---|---|---|",
    ]
    for row in rows:

        def cell(s):
            return s.replace("|", "\\|").replace("\n", " ")

        table.append(
            "| "
            + " | ".join(
                [
                    f"[{row['skill']}](../../{row['source']})",
                    cell(
                        row["classification"]
                        + " / "
                        + ", ".join(row["recipe_integration"]["adapters"])
                    ),
                    cell(row["inputs_units_outputs"]),
                    cell(row["recipe_integration"]["remaining_blockers"]),
                ]
            )
            + " |"
        )
    (DEST / "MATRIX.md").write_text("\n".join(table) + "\n")
    print(
        json.dumps(
            dict(
                skills=len(rows),
                scripts=sum(len(r["scripts"]) for r in rows),
                mcp_references=sum(len(r["mcp_entrypoints"]) for r in rows),
                unresolved=sum(
                    e["availability"] == "unresolved-reference"
                    for r in rows
                    for e in r["mcp_entrypoints"]
                ),
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
