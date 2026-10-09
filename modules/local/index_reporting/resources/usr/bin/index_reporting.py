#!/usr/bin/env python3
"""Format an LLM web service response as a MultiQC custom content HTML section.

Takes the raw JSON response saved by query_llm_service.py (OpenAI-style
chat completion or a custom {"result": "..."} response), extracts the text,
and renders it as an HTML table (when the LLM returned structured JSON) or
as preformatted text. An optional input JSON with structured candidates
(e.g. from prepare_demux_stats_for_llm) can be used to build richer tables.
"""

import argparse
import html
import json
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an LLM JSON response to a MultiQC custom content HTML."
    )
    parser.add_argument(
        "--response-json",
        "-r",
        required=True,
        help="Raw JSON response saved by query_llm_service.py.",
    )
    parser.add_argument(
        "--input-json",
        "-i",
        default="",
        help="Optional structured input JSON (e.g. CandidateRecoveries) used to build tables.",
    )
    parser.add_argument(
        "--output-mqc",
        "-o",
        required=True,
        help="Output HTML file path for the MultiQC custom content section.",
    )
    parser.add_argument(
        "--sample-name",
        default="LLM evaluation",
        help="Sample or run name displayed in the section (default: 'LLM evaluation').",
    )
    parser.add_argument(
        "--section-id",
        default="llm_evaluation",
        help="MultiQC section id (default: llm_evaluation).",
    )
    parser.add_argument(
        "--section-name",
        default="LLM evaluation",
        help="MultiQC section name (default: 'LLM evaluation').",
    )
    parser.add_argument(
        "--description",
        default="Automated LLM evaluation.",
        help="MultiQC section description.",
    )
    parser.add_argument(
        "--ok-message",
        default="No issues detected.",
        help="Message displayed when the evaluation reports no issues.",
    )
    return parser.parse_args()


def extract_content(response_data: dict) -> str:
    """Extract message content from OpenAI or alternative API response structure."""
    if "result" in response_data and isinstance(response_data["result"], str):
        return response_data["result"]
    if "choices" in response_data and len(response_data["choices"]) > 0:
        choice = response_data["choices"][0]
        if "message" in choice and "content" in choice["message"]:
            return choice["message"]["content"]
        if "text" in choice:
            return choice["text"]
    if "response" in response_data:
        return response_data["response"]
    if "content" in response_data:
        return response_data["content"]
    # Fallback to pretty-printed json if unknown format
    return json.dumps(response_data, indent=2)


def _format_cell(val, col_name: str = "", is_first: bool = False) -> str:
    col_lower = col_name.lower()
    if isinstance(val, list):
        if not val:
            return '<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; color:#64748b;">-</td>'
        count = len(val)
        items_str = ", ".join(str(x) for x in val)
        escaped_items = html.escape(items_str)
        if count <= 3:
            return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#334155;">{escaped_items}</td>'
        else:
            cell_content = (
                f'<div style="margin-bottom:4px;">'
                f'<span style="padding:2px 7px; border-radius:10px; background:#0284c7; color:#fff; font-size:11px; font-weight:600;">{count} items</span>'
                f'</div>'
                f'<div style="max-height:95px; overflow-y:auto; font-size:12px; line-height:1.45; color:#475569; word-break:break-word; background:#f8fafc; padding:6px 8px; border-radius:4px; border:1px solid #e2e8f0; font-family:monospace;">{escaped_items}</div>'
            )
            return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; min-width:220px;">{cell_content}</td>'

    if isinstance(val, str) and (
        val.startswith("<div")
        or val.startswith("<table")
        or val.startswith("<strong")
        or "samples affected" in val
        or "<span" in val
    ):
        min_w = " min-width:480px;" if ("table" in val or "samples affected" in val) else ""
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0;{min_w}">{val}</td>'

    val_str = str(val) if val is not None else ""
    escaped = html.escape(val_str)

    if "sample" in col_lower and "affected" not in col_lower and "barcode" not in col_lower:
        nowrap = " white-space:nowrap;" if len(val_str) < 30 else ""
        return f'<td style="padding:9px 12px; font-weight:600; color:#1e293b; border-bottom:1px solid #e2e8f0;{nowrap}">{escaped}</td>'
    elif "expected" in col_lower:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#475569; white-space:nowrap;">{escaped}</td>'
    elif any(w in col_lower for w in ["unassigned", "observed", "barcode"]) and "affected" not in col_lower:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; font-family:monospace; font-size:12.5px; color:#0284c7; font-weight:600; white-space:nowrap;">{escaped}</td>'
    elif "relation" in col_lower or "similarity" in col_lower or "issue" in col_lower:
        return f'<td style="padding:9px 12px; color:#b91c1c; font-weight:500; border-bottom:1px solid #e2e8f0;">{escaped}</td>'
    elif any(w in col_lower for w in ["recover", "read"]):
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0; white-space:nowrap;"><span style="display:inline-block; font-weight:600; color:#15803d; background:#ecfdf5; padding:3px 9px; border-radius:4px; border:1px solid #a7f3d0; font-size:12.5px;">{escaped}</span></td>'
    elif any(w in col_lower for w in ["action", "solution", "recommendation", "sheet"]):
        return f'<td style="padding:9px 12px; color:#0f766e; background:#f0fdfa; border-bottom:1px solid #e2e8f0;"><code style="color:#0f766e; background:none;">{escaped}</code></td>'
    else:
        return f'<td style="padding:9px 12px; border-bottom:1px solid #e2e8f0;">{escaped}</td>'


def _build_html_table(rows: list, headers: list) -> str:
    th_html = "".join(
        f'<th style="padding:10px 12px; font-weight:600; color:#334155; border-bottom:2px solid #cbd5e1;">{html.escape(h)}</th>'
        for h in headers
    )
    tbody_rows = []
    for r in rows:
        tds = []
        for i, val in enumerate(r):
            h_name = headers[i] if i < len(headers) else ""
            tds.append(_format_cell(val, col_name=h_name, is_first=(i == 0)))
        tbody_rows.append(f'<tr>{"".join(tds)}</tr>')
    tbody_html = "".join(tbody_rows)
    return f"""<div style="overflow-x:auto; margin-top:10px; margin-bottom:12px;">
  <table class="table table-bordered table-striped" style="width:100%; border-collapse:collapse; background:#ffffff; border:1px solid #e2e8f0; border-radius:6px; font-size:13.5px; text-align:left;">
    <thead><tr style="background:#f8fafc;">{th_html}</tr></thead>
    <tbody>{tbody_html}</tbody>
  </table>
</div>"""


def _extract_barcode_metrics(input_data: dict) -> tuple[dict, dict]:
    """Extract (undet_per_lane, barcode_reads) from Illumina Stats.json or similar."""
    if not input_data or not isinstance(input_data, dict):
        return {}, {}

    undet_per_lane = {}
    for cr in input_data.get("ConversionResults", []):
        l = cr.get("LaneNumber")
        undet = cr.get("Undetermined", {}).get("NumberReads", 0)
        if l is not None:
            undet_per_lane[l] = undet
            undet_per_lane[f"Lane {l}"] = undet

    barcode_reads = {}
    for ub in input_data.get("UnknownBarcodes", []):
        l = ub.get("Lane")
        for bc, count in ub.get("Barcodes", {}).items():
            if l is not None:
                barcode_reads[(l, bc)] = count
                barcode_reads[(f"Lane {l}", bc)] = count
            barcode_reads[bc] = count

    return undet_per_lane, barcode_reads


def _build_subtable_html(mappings: list, barcode_reads: dict = None) -> str:
    """Build an HTML sub-table for sample index mappings (Sample, Expected, Probable, Reads)."""
    if not mappings:
        return "-"

    barcode_reads = barcode_reads or {}
    tbody_lines = []
    has_reads = False
    mapping_regex = re.compile(
        r"^([^:]+):\s*([A-Za-z0-9_\-+]+)\s*(?:->|→|=>)\s*([A-Za-z0-9_\-+]+)(?:\s*\(([^)]+)\))?"
    )

    for idx, m in enumerate(mappings):
        bg_col = "#ffffff" if idx % 2 == 0 else "#f8fafc"
        if isinstance(m, dict):
            sid = m.get("SampleId", m.get("Sample", "-"))
            lanes = m.get("Lanes", m.get("Lane", ""))
            s_label = f"{sid} ({lanes})" if lanes and str(lanes) not in str(sid) else str(sid)
            frm = m.get("From", m.get("ExpectedIndex", "-"))
            to = m.get("To", m.get("ObservedBarcode", "-"))

            lane_num = None
            if lanes:
                match = re.search(r"\d+", str(lanes))
                if match:
                    lane_num = int(match.group(0))

            rds = None
            if "Reads" in m and str(m["Reads"]).isdigit():
                rds = int(m["Reads"])
            elif (lane_num, to) in barcode_reads:
                rds = barcode_reads[(lane_num, to)]
            elif to in barcode_reads:
                rds = barcode_reads[to]

            if rds is not None:
                has_reads = True
                rds_td = f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#15803d; font-weight:600; text-align:right; border-bottom:1px solid #f1f5f9; white-space:nowrap;">{rds:,}</td>'
            else:
                rds_td = '<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#94a3b8; text-align:center; border-bottom:1px solid #f1f5f9;">-</td>'

            tbody_lines.append(
                f'<tr style="background:{bg_col}; border-bottom:1px solid #f1f5f9;">'
                f'<td style="padding:5px 8px; font-weight:600; font-size:12px; color:#1e293b; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(s_label))}</td>'
                f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#475569; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(frm))}</td>'
                f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#0284c7; font-weight:600; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(to))}</td>'
                f'{rds_td}'
                f'</tr>'
            )
        else:
            match = mapping_regex.match(str(m).strip())
            if match:
                sid, exp_idx, obs_bc, rds = match.groups()
                lane_match = re.search(r"Lane\s*(\d+)", sid)
                lane_num = int(lane_match.group(1)) if lane_match else None
                rds_val = None
                if rds:
                    rds_clean = rds.replace("reads", "").replace("read", "").strip()
                    if rds_clean.isdigit():
                        rds_val = int(rds_clean)
                if rds_val is None:
                    if (lane_num, obs_bc) in barcode_reads:
                        rds_val = barcode_reads[(lane_num, obs_bc)]
                    elif obs_bc in barcode_reads:
                        rds_val = barcode_reads[obs_bc]

                if rds_val is not None:
                    has_reads = True
                    rds_td = f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#15803d; font-weight:600; text-align:right; border-bottom:1px solid #f1f5f9; white-space:nowrap;">{rds_val:,}</td>'
                else:
                    rds_td = '<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#94a3b8; text-align:center; border-bottom:1px solid #f1f5f9;">-</td>'

                tbody_lines.append(
                    f'<tr style="background:{bg_col}; border-bottom:1px solid #f1f5f9;">'
                    f'<td style="padding:5px 8px; font-weight:600; font-size:12px; color:#1e293b; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(sid))}</td>'
                    f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#475569; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(exp_idx))}</td>'
                    f'<td style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; color:#0284c7; font-weight:600; white-space:nowrap; border-bottom:1px solid #f1f5f9;">{html.escape(str(obs_bc))}</td>'
                    f'{rds_td}'
                    f'</tr>'
                )
            else:
                tbody_lines.append(
                    f'<tr style="background:{bg_col}; border-bottom:1px solid #f1f5f9;">'
                    f'<td colspan="4" style="padding:5px 8px; font-family:Menlo,Monaco,Consolas,monospace; font-size:12px; border-bottom:1px solid #f1f5f9;">{html.escape(str(m))}</td>'
                    f'</tr>'
                )

    reads_th = '<th style="padding:6px 8px; font-weight:600; font-size:11.5px; text-align:right; border-bottom:1px solid #cbd5e1;">Reads</th>' if has_reads else '<th style="padding:6px 8px; font-weight:600; font-size:11.5px; text-align:center; border-bottom:1px solid #cbd5e1;">Reads</th>'

    count = len(mappings)
    badge_html = (
        f'<div style="margin-bottom:6px;">'
        f'<span style="padding:3px 9px; border-radius:10px; background:#0284c7; color:#fff; font-size:11px; font-weight:600;">{count} samples affected</span>'
        f'</div>'
    )

    return f"""{badge_html}<div style="max-height:240px; overflow-y:auto; border:1px solid #cbd5e1; border-radius:4px; background:#ffffff;">
  <table style="width:100%; border-collapse:collapse; text-align:left;">
    <thead>
      <tr style="background:#f8fafc; color:#334155; position:sticky; top:0; z-index:1; border-bottom:2px solid #cbd5e1;">
        <th style="padding:6px 8px; font-weight:600; font-size:11.5px; border-bottom:1px solid #cbd5e1;">Sample</th>
        <th style="padding:6px 8px; font-weight:600; font-size:11.5px; border-bottom:1px solid #cbd5e1;">Expected Index</th>
        <th style="padding:6px 8px; font-weight:600; font-size:11.5px; border-bottom:1px solid #cbd5e1;">Probable Index</th>
        {reads_th}
      </tr>
    </thead>
    <tbody>
      {''.join(tbody_lines)}
    </tbody>
  </table>
</div>"""


def _parse_llm_json_to_table(parsed, input_data: dict = None) -> tuple[list, list] | None:
    if isinstance(parsed, dict):
        for wrapper in ["problematic_lanes", "lanes", "results", "issues", "diagnoses", "items", "recoveries"]:
            if wrapper in parsed and isinstance(parsed[wrapper], (list, dict)):
                parsed = parsed[wrapper]
                break

    rows = []
    if isinstance(parsed, dict):
        for lane, details in parsed.items():
            if isinstance(details, dict):
                issue = ""
                action = ""
                other_details = []
                for k, v in details.items():
                    k_lower = k.lower()
                    if any(w in k_lower for w in ["issue", "problem", "hypothesis", "diagnosis"]):
                        issue = str(v)
                    elif any(w in k_lower for w in ["action", "solution", "recommendation", "sheet"]):
                        action = str(v)
                    else:
                        other_details.append(f"{k}: {v}")
                if not issue and details:
                    first_k, first_v = next(iter(details.items()))
                    issue = f"{first_k}: {first_v}"
                if other_details:
                    action = f"{action} ({'; '.join(other_details)})" if action else "; ".join(other_details)
                rows.append([str(lane), issue or "-", action or "-"])
            elif isinstance(details, (str, int, float)):
                rows.append([str(lane), str(details), "-"])
        if rows:
            return ["Lane / Item", "Identified Issue", "Recommended Action"], rows

    elif isinstance(parsed, list):
        undet_per_lane, barcode_reads = _extract_barcode_metrics(input_data)

        # 1. Grouped issues (from LLM) with mappings or sample lists
        has_mappings_or_samples = any(
            isinstance(item, dict) and any(
                any(w in k.lower() for w in ["mapping", "sample", "barcode"])
                and isinstance(v, list)
                and len(v) > 0
                for k, v in item.items()
            )
            for item in parsed
        )

        if has_mappings_or_samples:
            grouped_rows = []
            for item in parsed:
                if not isinstance(item, dict):
                    continue

                issue = item.get("Issue", item.get("Diagnosis", item.get("Problem", "-")))
                rel = item.get("Relationship", "")
                gen_reads = item.get("RecoverableReads", item.get("PotentialRecoverableReads", item.get("Reads", "-")))
                action = item.get("Action", item.get("RecommendedAction", "-"))

                # Combine Issue and Relationship nicely
                if rel and issue and rel.strip().lower() != issue.strip().lower() and issue != "-":
                    issue_label = f"<strong>{html.escape(str(issue))}</strong><br><span style='color:#64748b; font-size:12px;'>Relationship: <span style='color:#b91c1c; font-weight:600;'>{html.escape(str(rel))}</span></span>"
                elif rel:
                    issue_label = f"<strong>{html.escape(str(rel))}</strong>"
                elif issue:
                    issue_label = f"<strong>{html.escape(str(issue))}</strong>"
                else:
                    issue_label = "-"

                # Extract mappings or sample list
                mappings = []
                for k, v in item.items():
                    if "mapping" in k.lower() and isinstance(v, list) and v:
                        mappings = v
                        break
                if not mappings:
                    for k, v in item.items():
                        if any(w in k.lower() for w in ["sample", "barcode"]) and isinstance(v, list) and v:
                            mappings = v
                            break

                # Subtable HTML with Sample, Indice Atteso, Indice Probabile, Reads
                subtable_html = _build_subtable_html(mappings, barcode_reads=barcode_reads)

                # Recoverable reads display
                tot_undet = sum(v for k, v in undet_per_lane.items() if isinstance(k, int))
                if str(gen_reads).isdigit():
                    rds_val = int(gen_reads)
                    if tot_undet > 0:
                        pct = (rds_val / tot_undet) * 100
                        rds_str = f"{rds_val:,} ({pct:.2f}%)"
                    else:
                        rds_str = f"{rds_val:,}"
                elif gen_reads and gen_reads != "-":
                    rds_str = str(gen_reads)
                else:
                    rds_str = "-"

                grouped_rows.append([issue_label, rds_str, str(action), subtable_html])

            if grouped_rows:
                headers = [
                    "Identified Issue",
                    "Recoverable Reads",
                    "Recommended Action",
                    "Affected Samples",
                ]
                return headers, grouped_rows

        # 2. Check if list items represent direct mapping objects
        is_direct_mapping = any(
            isinstance(item, dict) and any(
                any(w in k.lower() for w in ["expectedindex", "observedbarcode", "unknownbarcode", "candidatesample"])
                for k in item.keys()
            )
            for item in parsed
        )
        if is_direct_mapping:
            unpacked_rows = []
            for item in parsed:
                if isinstance(item, dict):
                    sid = item.get("Sample", item.get("CandidateSample", item.get("SampleId", "-")))
                    exp_idx = item.get("ExpectedIndex", "-")
                    obs_bc = item.get("ObservedBarcode", item.get("UnknownBarcode", item.get("FoundInUnassigned", "-")))
                    rel = item.get("Relationship", item.get("Issue", "-"))
                    rds = item.get("RecoverableReads", item.get("Reads", "-"))
                    action = item.get("Action", item.get("RecommendedAction", "-"))
                    rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
                    unpacked_rows.append([str(sid), str(exp_idx), str(obs_bc), str(rel), rds_str, str(action)])
            if unpacked_rows:
                headers = [
                    "Candidate Sample",
                    "Expected Index",
                    "Found in Unassigned",
                    "Relationship / Similarity",
                    "Recoverable Reads",
                    "Recommended Action",
                ]
                return headers, unpacked_rows

        # 3. Fallback for simple lane or issue items
        has_recovery = any(
            isinstance(item, dict) and any(
                any(w in k.lower() for w in ["recover", "read_count", "potential_reads", "reads"])
                for k in item.keys()
            )
            for item in parsed
        )
        for item in parsed:
            if isinstance(item, dict):
                lane = ""
                issue = ""
                action = ""
                recovery = ""
                for k, v in item.items():
                    k_lower = k.lower()
                    if any(w in k_lower for w in ["lane", "sample", "id"]):
                        lane = str(v) if not isinstance(v, list) else ", ".join(str(x) for x in v)
                    elif any(w in k_lower for w in ["issue", "problem", "hypothesis", "diagnosis"]):
                        issue = str(v)
                    elif any(w in k_lower for w in ["action", "solution", "recommendation", "sheet"]):
                        action = str(v)
                    elif any(w in k_lower for w in ["recover", "read_count", "potential_reads", "reads"]):
                        recovery = str(v)
                if not lane and item:
                    lane = f"Lane {len(rows)+1}"
                if has_recovery:
                    rows.append([lane or "-", issue or "-", recovery or "-", action or "-"])
                else:
                    rows.append([lane or "-", issue or "-", action or "-"])
            elif isinstance(item, str):
                rows.append([f"Item {len(rows)+1}", str(item), "-"])
        if rows:
            headers = (
                ["Lane / Item", "Identified Issue", "Potential Recoverable Reads", "Recommended Action"]
                if has_recovery
                else ["Lane / Item", "Identified Issue", "Recommended Action"]
            )
            return headers, rows

    return None


def generate_multiqc_html(
    content: str,
    output_path: str,
    sample_name: str = "LLM evaluation",
    input_data: dict = None,
    section_id: str = "llm_evaluation",
    section_name: str = "LLM evaluation",
    description: str = "Automated LLM evaluation.",
    ok_message: str = "No issues detected.",
):
    """Generate a MultiQC custom content HTML file."""
    clean_text = content.strip()
    if clean_text.startswith("```"):
        lines = clean_text.splitlines()
        if lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].startswith("```"):
            lines = lines[:-1]
        clean_text = "\n".join(lines).strip()

    # Parse LLM JSON if possible
    parsed = None
    try:
        parsed = json.loads(clean_text)
    except Exception:
        pass

    # Check if structured data is provided via input_data
    has_structured_data = bool(
        input_data
        and (
            input_data.get("CandidateRecoveries")
            or input_data.get("UnmatchedZeroReadSamples")
            or input_data.get("UnassignedAbundantBarcodes")
        )
    )

    if has_structured_data:
        headers = [
            "Candidate Sample",
            "Expected Index",
            "Found in Unassigned",
            "Relationship / Similarity",
            "Recoverable Reads",
            "Recommended Action",
        ]
        rows = []

        llm_unmatched_action = ""
        llm_unassigned_action = ""
        if parsed and isinstance(parsed, list):
            for item in parsed:
                if isinstance(item, dict):
                    iss = str(item.get("Issue", "")).lower()
                    act = str(item.get("Action", item.get("RecommendedAction", "")))
                    if any(w in iss for w in ["zero", "failure", "massive", "mismatch"]):
                        llm_unmatched_action = act
                    elif any(w in iss for w in ["unassigned", "unknown", "spike-in", "omitted"]):
                        llm_unassigned_action = act

        # 1. Candidate Recoveries
        for c in input_data.get("CandidateRecoveries", []):
            sid = c.get("CandidateSample", "-")
            lane = c.get("Lane", "")
            s_label = f"{sid} ({lane})" if lane else sid
            exp_idx = c.get("ExpectedIndex", "-")
            obs_bc = c.get("UnknownBarcode", "-")
            rel = c.get("Relationship", "-")
            rds = c.get("Reads", 0)
            rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
            action = c.get("RecommendedAction", "-")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        # 2. Unmatched Zero-Read Samples
        for s in input_data.get("UnmatchedZeroReadSamples", []):
            sid = s.get("SampleId", "-")
            lane = s.get("Lane", "")
            s_label = f"{sid} ({lane})" if lane else sid
            exp_idx = s.get("ExpectedIndex", "-")
            obs_bc = s.get("FoundInUnassigned", "None found")
            rel = s.get("Relationship", "Zero reads (no match in undetermined)")
            rds = s.get("Reads", 0)
            rds_str = str(rds)
            action = llm_unmatched_action or s.get("RecommendedAction", "Verify sample sheet index against library preparation design")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        # 3. Unassigned Abundant Barcodes
        for u in input_data.get("UnassignedAbundantBarcodes", []):
            lane = u.get("Lane", "")
            s_label = f"[Unassigned Barcode] ({lane})" if lane else "[Unassigned Barcode]"
            exp_idx = "-"
            obs_bc = u.get("UnknownBarcode", "-")
            rel = u.get("Relationship", "Abundant unassigned barcode")
            rds = u.get("Reads", 0)
            rds_str = f"{int(rds):,}" if str(rds).isdigit() else str(rds)
            action = llm_unassigned_action or u.get("RecommendedAction", "Check for omitted sample or spike-in in library")
            rows.append([s_label, exp_idx, obs_bc, rel, rds_str, action])

        is_ok = len(rows) == 0
        badge_color = "#28a745" if is_ok else "#dc3545"
        status_label = "PASS" if is_ok else "ATTENTION NEEDED"

        if is_ok:
            body_display = f"""<div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #22c55e; border-radius:4px; padding:12px 16px; margin:10px 0; color:#15803d; font-size:13.5px;">
  <strong>✓ {html.escape(section_name)}:</strong> {html.escape(ok_message)}
</div>"""
        else:
            summary_boxes_html = ""
            if parsed and isinstance(parsed, list):
                boxes = []
                for item in parsed:
                    if isinstance(item, dict):
                        iss = item.get("Issue", "")
                        act = item.get("Action", item.get("RecommendedAction", ""))
                        if iss or act:
                            boxes.append(f"""<div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:4px; padding:10px 14px; margin:8px 0; font-size:13px; color:#991b1b;">
  <strong>Identified Issue:</strong> {html.escape(str(iss))}<br>
  <strong>Recommendation:</strong> {html.escape(str(act))}
</div>""")
                if boxes:
                    summary_boxes_html = "\n".join(boxes)

            table_html = _build_html_table(rows, headers)
            body_display = f"{summary_boxes_html}\n{table_html}" if summary_boxes_html else table_html

    else:
        # Fallback when structured data is not available
        is_ok = (
            "no problems" in clean_text.lower()
            or "no issues" in clean_text.lower()
            or "all correct" in clean_text.lower()
        )
        badge_color = (
            "#28a745"
            if is_ok
            else "#e0a800"
            if "warn" in clean_text.lower()
            else "#dc3545"
        )
        status_label = "PASS" if is_ok else "ATTENTION NEEDED"
        escaped_content = html.escape(clean_text)

        if is_ok:
            body_display = f"""<div style="background:#f0fdf4; border:1px solid #bbf7d0; border-left:4px solid #22c55e; border-radius:4px; padding:12px 16px; margin:10px 0; color:#15803d; font-size:13.5px;">
  <strong>✓ {html.escape(section_name)}:</strong> {html.escape(ok_message)}
</div>"""
        elif parsed is not None:
            table_data = _parse_llm_json_to_table(parsed, input_data=input_data)
            pretty_json = html.escape(json.dumps(parsed, indent=2))

            summary_boxes_html = ""
            if parsed and isinstance(parsed, list):
                boxes = []
                for item in parsed:
                    if isinstance(item, dict):
                        iss = item.get("Issue", "")
                        act = item.get("Action", item.get("RecommendedAction", ""))
                        if iss or act:
                            boxes.append(f"""<div style="background:#fef2f2; border:1px solid #fecaca; border-left:4px solid #ef4444; border-radius:4px; padding:10px 14px; margin:8px 0; font-size:13px; color:#991b1b;">
  <strong>Identified Issue:</strong> {html.escape(str(iss))}<br>
  <strong>Recommendation:</strong> {html.escape(str(act))}
</div>""")
                if boxes:
                    summary_boxes_html = "\n".join(boxes)

            if table_data:
                headers, rows = table_data
                table_html = _build_html_table(rows, headers)
                body_display = f"{summary_boxes_html}\n{table_html}" if summary_boxes_html else table_html
            else:
                body_display = f"{summary_boxes_html}\n<pre style=\"background:#f8f9fa; border:1px solid #e9ecef; border-radius:4px; padding:12px; margin-top:8px; font-size:13px; font-family:Menlo,Monaco,Consolas,monospace;\"><code>{pretty_json}</code></pre>"
        else:
            body_display = f'<pre style="white-space:pre-wrap; background:#f8f9fa; border:1px solid #e9ecef; border-radius:4px; padding:12px; margin-top:8px; font-size:13px; font-family:Menlo,Monaco,Consolas,monospace;">{escaped_content}</pre>'

    html_content = f"""<!--
id: '{section_id}'
section_name: '{section_name}'
description: '{description}'
-->
<div style="font-family: inherit; margin: 15px 0;">
  <div style="display: flex; align-items: center; margin-bottom: 8px;">
    <h4 style="margin: 0; font-size: 1.05em; font-weight: 600;">{html.escape(sample_name)}</h4>
    <span style="margin-left: 12px; padding: 2px 8px; border-radius: 4px; font-size: 0.78em; font-weight: bold; color: #fff; background-color: {badge_color};">
      {status_label}
    </span>
  </div>
  {body_display}
</div>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)


def main():
    args = parse_args()

    try:
        with open(args.response_json, "r", encoding="utf-8") as f:
            response_json = json.load(f)
    except Exception as e:
        print(f"Error reading response JSON '{args.response_json}': {e}", file=sys.stderr)
        sys.exit(1)

    input_data = None
    if args.input_json:
        try:
            with open(args.input_json, "r", encoding="utf-8") as f:
                input_data = json.load(f)
        except Exception as e:
            print(f"Warning: cannot read input JSON '{args.input_json}': {e}", file=sys.stderr)

    content = extract_content(response_json)
    generate_multiqc_html(
        content,
        args.output_mqc,
        sample_name=args.sample_name,
        input_data=input_data,
        section_id=args.section_id,
        section_name=args.section_name,
        description=args.description,
        ok_message=args.ok_message,
    )
    print(f"MultiQC section saved to: {args.output_mqc}")


if __name__ == "__main__":
    main()
